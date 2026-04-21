#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include "lbm.h"
#define Q 9
#define save_iter 50
#define disp_iter 10
#define N 2

// Ancillary functions
void ade_solver(double dt, double p_fluc, 
                const double lambdas[N], const double alpha[N], const double beta[N],
                const double phi_prev[N], const double psi_prev[N], const double chi_prev[N],
                double phi_out[N], double psi_out[N], double chi_out[N]) {
    
    for (int i = 0; i < N; i++) {
        phi_out[i] = phi_prev[i] + dt * (p_fluc + lambdas[i] * phi_prev[i]);
        psi_out[i] = psi_prev[i] + dt * (p_fluc + alpha[i] * psi_prev[i] - beta[i] * chi_prev[i]);
        chi_out[i] = chi_prev[i] + dt * (beta[i] * psi_prev[i] + alpha[i] * chi_prev[i]);
    }
}

double compute_acoustic_velocity(double p_fluc, double dt, double Zinf,
                                 const double A[N], const double B[N], const double C[N],
                                 const double lambdas[N], const double alpha[N], const double beta[N],
                                 double phi[N], double psi[N], double chi[N]) {

    // Update states
    ade_solver(dt, p_fluc, lambdas, alpha, beta, phi, psi, chi, phi, psi, chi);

    // Compute wall-normal velocity
    double dot_A_phi = 0.0;
    double dot_B_psi = 0.0;
    double dot_C_chi = 0.0;

    for (int i = 0; i < N; i++) {
        dot_A_phi += A[i] * phi[i];
        dot_B_psi += B[i] * psi[i];
        dot_C_chi += C[i] * chi[i];
    }

    double v = Zinf * p_fluc - dot_A_phi - 2.0 * (dot_B_psi + dot_C_chi);

    return v;
}

int load_impedance_data(const char *filename, double *zinf_ptr, double a[], 
                    double l[], double b[], double c[], 
                    double al[], double be[]) {
    
    FILE *fptr = fopen(filename, "r");

    // Error opening file
    if (fptr == NULL) {
        return 1; 
    }

    // Read single value into the memory address provided
    if (fscanf(fptr, "%lf", zinf_ptr) != 1) return 1;

    // Read array values
    for (int i = 0; i < N; i++) if (fscanf(fptr, "%lf", &a[i]) != 1) return 1;
    for (int i = 0; i < N; i++) if (fscanf(fptr, "%lf", &l[i]) != 1) return 1;
    for (int i = 0; i < N; i++) if (fscanf(fptr, "%lf", &b[i]) != 1) return 1;
    for (int i = 0; i < N; i++) if (fscanf(fptr, "%lf", &c[i]) != 1) return 1;
    for (int i = 0; i < N; i++) if (fscanf(fptr, "%lf", &al[i]) != 1) return 1;
    for (int i = 0; i < N; i++) if (fscanf(fptr, "%lf", &be[i]) != 1) return 1;

    fclose(fptr);

    return 0;
}

// Main function
int main(int argc, char *argv[]){

    // Check if the correct number of arguments is provided
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <p_a> <f_exc>\n", argv[0]);
        fprintf(stderr, "Example: %s 89 800\n", argv[0]);
        return 1;
    }

    // UI parameters (pressure amplitude [Pa] and excitation frequency [Hz])
    double p_a = atof(argv[1]);                   // Pressure amplitude [Pa] (89 -> 130 dB, 503 -> 145 dB)
    double f_exc = atof(argv[2]);                 // Excitation frequency [Hz] (800, 1000, 1400, 2000)

    // Computational domain (physical units)
    int N_o = 2;                                  // Number of orifices [-]
    double Ly = 0.009906;                         // Domain length along y [m]
    double y1 = 0.002528;                         // Orifice 1 location along y [m]
    double y2 = 0.007478;                         // Orifice 2 location along y [m]
    double d = 0.00117;                           // Orifice diameter [m]
    int ac_cycles = 10;                           // Number of acoustic cycles [-]
    int damping_cycles = 2;                       // Number of damping cycles [-]
    int start_cycles = 1;                         // Number of starting cycles [-]
    double c0 = 340.0;                            // Sound speed [m/s]
    double tau = 0.5005;                          // Relaxation time [-]
    double nu = 1.5e-5;                           // Kinematic viscosity [m^2/s]
    double rho_0 = 1.184;                         // Uniform density at rest [kg/m^3]

    // Conversion factors (physical to lattice units)
    double cs = 1.0 / sqrt(3.0);                  // Speed of sound [m/s]
    double dx = 3 * nu * cs / (c0 * (tau - 0.5)); // Voxel size [m]
    double crho = 1.0;                            // Lattice density conversion factor [kg/m^3]
    double dt = dx * dx * (tau - 0.5) / (3 * nu); // Time step [s] dx * dx * (tau - 0.5) / (3 * nu)
    double cu = dx / dt;                          // Lattice velocity conversion factor [m/s]
    double cp = crho * cu * cu;                   // Lattice pressure conversion factor [Pa]

    // Forcing (lattice units)
    double Fx = 0.0;                              // Body force x component [m / s^2]
    double Fy = 0.0;                              // Body force y component [m / s^2]

    // Initial conditions (lattice units)
    int init_iter_max = 100;                      // Maximum number of iterations in Mei's algorithm [-]
    double tol_rho = 1e-10;                       // Tolerance for density convergence in Mei's algorithm [-]

    // Compute derivate quantities
    int tot_cycles = start_cycles + ac_cycles + damping_cycles;
    double lambda = c0 / f_exc;                   // Wavelength along x [m]
    double Lx = tot_cycles * lambda;              // Domain length along x [m]
    double kx = 2 * M_PI / lambda;                // Wavenumber along x [1/m]
    double rho_a = (p_a / cp) / (cs * cs);        // Density amplitude [kg/m^3]
    double x_end = Lx - start_cycles*lambda;      // Wave packet end [m]
    double x_start = x_end - ac_cycles*lambda;    // Wave packet start [m]
    int Nx = ceil(Lx / dx);                       // Number of points along x [-]
    int Ny = ceil(Ly / dx);                       // Number of points along y [-]
    int Nt = ceil(tot_cycles / (f_exc*dt));       // Number of time steps [-]

    // Read impedance fitting data
    double Zinf;
    double A[N], lambdas[N], B[N], C[N], alpha[N], beta[N];
    char filename_impedance[256];

    snprintf(filename_impedance, sizeof(filename_impedance),"impedance_%gPa.txt", p_a);
    if (load_impedance_data(filename_impedance, &Zinf, A, lambdas, B, C, alpha, beta) == 0) {
        printf("Impedance data loaded successfully.\n");
    } else {
        printf("Failed to load impedance data.\n");
    }

    // Initialize ADE states
    double phi[N] = {0, 0};
    double psi[N] = {0, 0};
    double chi[N] = {0, 0};

    // Boundary conditions
    int use_IBB = 0;                              // Use IBB (1: yes, 0: no) [-]
    int num_boundaries = 3;                       // Number of boundary conditions [-]

    Boundary boundaries[num_boundaries];

    boundaries[0].apply = periodic_y;

    double r0[Ny];
    for (int j = 0; j < Ny; j++){
        r0[j] = rho_0;
    }

    boundaries[1].apply = pressure_inlet;
    boundaries[1].index = 0;
    boundaries[1].val1 = r0;

    double dmy[Ny];
    for (int j = 0; j < Ny; j++) {
        dmy[j] = 0.0;
    } 

    boundaries[2].apply = velocity_outlet_regularized;
    boundaries[2].index = Nx-1;
    boundaries[2].val1 = dmy;
    boundaries[2].val2 = dmy;

    // Check on max Mach number in lattice units
    double Ma_max = rho_a / rho_0;
    if (Ma_max > 0.2) {
        printf("Warning: Maximum Mach number in lattice units > 0.2 (%.2f)\n", Ma_max);
    }

    // Display simulation info
    printf("Lattice size: %g - Time step: %g - Maximum lattice Mach: %g\n", dx, dt, Ma_max);
    printf("Pressure amplitude: %g Pa - Excitation frequency: %g Hz\n", p_a, f_exc);

    // Check if there is a periodic BC along x and/or y
    int isperiodic_x = 0;
    int isperiodic_y = 0;
    for (int i = 0; i < num_boundaries; i++) {
        if (boundaries[i].apply == periodic_x) isperiodic_x = 1;
        if (boundaries[i].apply == periodic_y) isperiodic_y = 1;
    }

    // Lattice velocity x and y components
    int cx[Q] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
    int cy[Q] = {0, 0, 1, 0, -1, 1, 1, -1, -1};

    // Lattice weights (D2Q9)
    double w[Q] = {4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};

    // Indices of opposite directions for D2Q9 lattice
    int opp[Q] = {0, 3, 4, 1, 2, 7, 8, 5, 6};

    // Solid mask (staircase approximation) and signed distance (IBB)
    int (*solid_mask)[Nx] = malloc(Ny * sizeof *solid_mask);
    double (*phi_solid)[Nx] = malloc(Ny * sizeof *phi_solid);
    for (int j = 0; j < Ny; j++){
        for (int i = 0; i < Nx; i++){
            solid_mask[j][i] = 0;
            phi_solid[j][i] = 0.0;
        }
    }

    // BGK collision operator coefficient (Omega = omega*(f - f_eq))
    double (*omega_eff)[Nx] = malloc(Ny * sizeof *omega_eff);
    for (int j = 0; j < Ny; j++){
        for (int i = 0; i < Nx; i++){
            omega_eff[j][i] = 1.0 / tau;
        }
    }

    // Flow field and particle distribution function initialization
    double (*U_in)[Nx] = malloc(Ny * sizeof *U_in);
    double (*V_in)[Nx] = malloc(Ny * sizeof *V_in);
    double (*rho_in)[Nx] = malloc(Ny * sizeof *rho_in);
    
    double (*rho)[Nx] = malloc(Ny * sizeof *rho);
    double (*u)[Nx] = malloc(Ny * sizeof *u);
    double (*v)[Nx] = malloc(Ny * sizeof *v);

    double (*f)[Nx][Q] = malloc(Ny * sizeof *f);
    double (*f_new)[Nx][Q] = malloc(Ny * sizeof *f_new);

    for (int j = 0; j < Ny; j++){
        for (int i = 0; i < Nx; i++){ 
            double x = i * dx;

            double envelope = ((x >= x_start) & (x <= x_end)) ? 1.0 : 0.0;
            double carrier  = sin(kx * x);

            double rho_p = rho_a * envelope * carrier;
            rho_in[j][i] = rho_0 + rho_p;

            U_in[j][i] = cs * rho_p / rho_0;
            V_in[j][i] = 0.0;
        }
    }
    mei_initialization(Nx, Ny, Q, f, f_new, rho_in, U_in, V_in, rho, u, v, cx, cy, w, omega_eff, isperiodic_x, isperiodic_y, init_iter_max, tol_rho, solid_mask);

    // Define and initialize the lift and drag forces on the surface (physical units)
    double L = 0.0;
    double D = 0.0;

    double p_fluc = 0.0;
    double v_ac = 0.0;
    double vBC[Ny];

    // Ensure that the "sol" directory exists and, if not, create it
    ensure_directory_exists("sol");

    // Open a file for writing
    char filename_tdibc[256];
    snprintf(filename_tdibc, sizeof(filename_tdibc),"sol/tdibc_%gHz_%gPa.txt", f_exc, p_a);
    FILE *fp = fopen(filename_tdibc, "w");

    // Main LBM loop
    for (int it = 0; it < Nt; it++){

        // Print the timestep
        if (it % disp_iter == 0) printf("Step %d of %d - p_fluc = %g - v_ac = %g\n ", it, Nt, p_fluc, v_ac);

        // Save acoustic pressure and velocity
        fprintf(fp, "%lf %lf %lf\n", it * dt, p_fluc, v_ac);

        // Compute surface averaged pressure fluctuations (physical units)
        double rho_surf = 0.0;
        for (int j = 0; j < Ny; j++){
            rho_surf += rho[j][Nx - 1];
        }
        rho_surf /= Ny;
        p_fluc = cs * cs * (rho_surf - rho_0) * cp;

        // Compute acoustic velocity solving ADE IVP (physical units)
        v_ac = compute_acoustic_velocity(p_fluc, dt, Zinf, A, B, C, lambdas, alpha, beta, phi, psi, chi);

        double temporary_velocity = (4 * Ly * Ly) / (N_o * d * d) * (v_ac / cu);
        if ((temporary_velocity / cs) > 0.2) {
            printf("Warning, max Mach number greater than 0.2 (%g)!\n", temporary_velocity / cs);
        }

        // Fill the local values of the acoustic velocity on the cavity surface (lattice units)
        for (int j = 0; j < Ny; j++){
            double y = j * dx;

            if ((y >= y1 - d/2 && y <= y1 + d/2) | (y >= y2 - d/2 && y <= y2 + d/2)) {
                vBC[j] = temporary_velocity;
            } else {
                vBC[j] = 0.0;
            }
        }

        // Apply acoustic velocity TDIM BC (lattice units)
        boundaries[2].val1 = vBC;

        // Execute the streaming and colliding steps
        main_lbm(Nx,Ny,Q,f,f_new,rho,u,v,cx,cy,w,opp,omega_eff,solid_mask,boundaries,num_boundaries,isperiodic_x,isperiodic_y,Fx,Fy,&L,&D,phi_solid,use_IBB);

        // Swap f and f_new
        double (*temp_ptr)[Nx][Q] = f;
        f = f_new;
        f_new = temp_ptr;

        // Write solution (rho, u, v) on a ".vtk" file each save_iter iterations
        if (it % save_iter == 0){
            char filename[256];
            snprintf(filename, sizeof(filename),"sol/fields_%gHz_%gPa_%05d.vtk", f_exc, p_a, it);
            write_vtk_binary_2D(filename, Nx, Ny, dx, u, v, rho, cu, crho);
        }

    }

    fclose(fp);

    // Free the memory
    free(U_in);
    free(V_in);
    free(rho_in);
    free(solid_mask);
    free(phi_solid);
    free(omega_eff);
    free(rho);
    free(u);
    free(v);
    free(f);
    free(f_new);

    return 0;

}
