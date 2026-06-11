#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include "lbm.h"
#define Q 9
#define disp_iter 10

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
    double d = 0.00117;                           // Orifice diameter (CHAR LENGTH) [m]
    double Ly = 0.009906;                         // Domain length along y [m]
    double tc = 5.500e-4;                         // Facesheet thickness [m]
    double zt = 3.865e-2;                         // Cavity length [m]
    int res = 20;                                 // Resolution (voxels per char length) [-]
    int ac_cycles = 10;                           // Number of acoustic cycles [-]
    int damping_cycles = 2;                       // Number of damping cycles [-]
    int start_cycles = 0;                         // Number of starting cycles [-]
    double c0 = 340.0;                            // Sound speed [m/s]
    double tau_target = 0.5063;                   // Relaxation time [-]
    double nu = 1.5e-5;                           // Kinematic viscosity [m^2/s]
    double rho_0 = 1.184;                         // Uniform density at rest [kg/m^3]

    // Conversion factors (physical to lattice units)
    double crho = 1.0;                            // Lattice density conversion factor [kg/m^3]
    double cs = 1.0 / sqrt(3.0);                  // Speed of sound [m/s]
    double dx = d / res;                          // Voxel size [m]
    double tau = 0.5 + 3 * nu * cs / (c0 * dx);   // Relaxation time [-]
    double tau_turb = tau_target - tau;           // Turbulence relaxation time [-]
    double dt = dx * dx * (tau - 0.5) / (3 * nu); // Time step [s]
    double cu = dx / dt;                          // Lattice velocity conversion factor [m/s]
    double cp = crho * cu * cu;                   // Lattice pressure conversion factor [Pa]

    // Forcing (lattice units)
    double Fx = 0.0;                              // Body force x component [m / s^2]
    double Fy = 0.0;                              // Body force y component [m / s^2]

    // Initial conditions (lattice units)
    int init_iter_max = 10;                       // Maximum number of iterations in Mei's algorithm [-]
    double tol_rho = 1e-10;                       // Tolerance for density convergence in Mei's algorithm [-]

    // Compute derivate quantities
    int tot_cycles = start_cycles + ac_cycles + damping_cycles;
    int sim_cycles = start_cycles + ac_cycles; 
    double lambda = c0 / f_exc;                   // Wavelength along x [m]
    double Lx = tot_cycles * lambda + zt + 2 * tc;// Domain length along x [m]
    double kx = 2 * M_PI / lambda;                // Wavenumber along x [1/m]
    double rho_a = (p_a / cp) / (cs * cs);        // Density amplitude [kg/m^3]
    double x_end = tot_cycles * lambda - start_cycles*lambda;      // Wave packet end [m]
    double x_start = x_end - ac_cycles*lambda;    // Wave packet start [m]
    double dt_eff = 1.0 / (2 * f_exc);            // Effective time step [s]
    double y1 = 0.255198869 * Ly;                 // Orifice 1 location along y [m]
    double y2 = 0.754896023 * Ly;                 // Orifice 2 location along y [m]
    int Nx = ceil(Lx / dx);                       // Number of points along x [-]
    int Ny = ceil((Ly + 2 * tc) / dx);            // Number of points along y [-]
    int Nt = ceil(sim_cycles / (f_exc * dt));     // Number of time steps [-]
    int init_save = floor((start_cycles + 0.25) / (f_exc * dt));   // Number of iterations to start saving [-]
    int save_iter = floor(dt_eff / dt);           // Every how many iterations to save [-]

    // Boundary conditions
    int isperiodic_x = 0;
    int isperiodic_y = 1;

    double r0[Ny];
    for (int j = 0; j < Ny; j++){
        r0[j] = rho_0;
    }

    // Check on max Mach number in lattice units
    printf("\n");
    double Ma_max = rho_a / rho_0;
    if (Ma_max > 0.2) printf("Warning: Maximum Mach number in lattice units > 0.2 (%.2f)\n", Ma_max);

    // Check on points per wavelength
    double points_per_wavelength = lambda / dx;
    if (points_per_wavelength < 20) {
        printf("Warning: Less than 20 points per wavelength (%.2f)\n", points_per_wavelength);
    } else {
        printf("Nx: %d, Ny: %d, Points per wavelength: %.2f\n", Nx, Ny, points_per_wavelength);
    }

    // Display simulation info
    printf("Relaxation time: %g - Turbulence relaxation time: %g\n", tau, tau_turb);
    printf("Lattice size: %g - Time step: %g - Maximum lattice Mach: %g\n", dx, dt, Ma_max);
    printf("Pressure amplitude: %g Pa - Excitation frequency: %g Hz\n", p_a, f_exc);

    // Lattice velocity x and y components
    int cx[Q] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
    int cy[Q] = {0, 0, 1, 0, -1, 1, 1, -1, -1};

    // Lattice weights (D2Q9)
    double w[Q] = {4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};

    // Indices of opposite directions for D2Q9 lattice
    int opp[Q] = {0, 3, 4, 1, 2, 7, 8, 5, 6};

    // Solid mask (staircase approximation) and signed distance (IBB)
    int (*solid_mask)[Nx] = malloc(Ny * sizeof *solid_mask);
    double x_duct_end = tot_cycles * lambda;
    for (int j = 0; j < Ny; j++){
        double y = j * dx;
        for (int i = 0; i < Nx; i++){
            double x = i * dx;
            if (((x >= x_duct_end) & (x <= x_duct_end + tc)) | ((x >= Lx - tc) & (x <= Lx))) {
                solid_mask[j][i] = 1;
            } else if ((x >= x_duct_end) & (((y >= 0) & (y <= tc)) | ((y >= Ly + tc) & (y <= Ly + 2 * tc)))) {
                solid_mask[j][i] = 1;
            } else {
                solid_mask[j][i] = 0;
            }

            if (((x >= x_duct_end) & (x <= x_duct_end + tc)) & (((y >= tc + y1 - d/2) & (y <= tc + y1 + d/2)) | ((y >= tc + y2 - d/2) & (y <= tc + y2 + d/2)))) {
                solid_mask[j][i] = 0;
            }
        }
    }

    // BGK collision operator coefficient (Omega = omega*(f - f_eq))
    double (*omega_eff)[Nx] = malloc(Ny * sizeof *omega_eff);
    double tau_eff = tau + tau_turb;
    for (int j = 0; j < Ny; j++){
        for (int i = 0; i < Nx; i++){
            omega_eff[j][i] = 1.0 / tau_eff;
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

            double envelope = ((x >= x_start) && (x <= x_end)) ? 1.0 : 0.0;
            double carrier  = sin(kx * x);

            double rho_p = rho_a * envelope * carrier;
            rho_in[j][i] = rho_0 + rho_p;

            U_in[j][i] = cs * rho_p / rho_0;
            V_in[j][i] = 0.0;
        }
    }
    mei_initialization(Nx, Ny, Q, f, f_new, rho_in, U_in, V_in, rho, u, v, cx, cy, w, omega_eff, isperiodic_x, isperiodic_y, init_iter_max, tol_rho, solid_mask);

    // Ensure that the "sol" directory exists and, if not, create it
    ensure_directory_exists("sol");

    // Open a file for writing the tdibc (pressure and velocity on the cavity)
    char filename_tdibc[256];
    snprintf(filename_tdibc, sizeof(filename_tdibc),"sol/tdibc_%gHz_%gPa.txt", f_exc, p_a);
    FILE *fp = fopen(filename_tdibc, "w");

    // Write initial solution
    char filename[256];
    snprintf(filename, sizeof(filename),"sol/fields_%gHz_%gPa_init.vtk", f_exc, p_a);
    write_vtk_binary_2D(filename, Nx, Ny, dx, u, v, rho, cu, crho);

    double rho_surf = 0.0;
    double p_ac = 0.0;
    double v_ac = 0.0;
    int i_surf = floor(x_duct_end / dx) - 1;

    // Main LBM loop
    for (int it = 0; it < Nt; it++){

        // Print the timestep
        if (it % disp_iter == 0) printf("Step %d of %d - p_ac = %g - v_ac = %g\n ", it, Nt, p_ac, v_ac);

        // Save acoustic quantities on the liner surface
        fprintf(fp, "%lf %lf %lf\n", it * dt, p_ac, v_ac);

        // Collision
        collision(Nx,Ny,Q,f,rho,u,v,omega_eff,solid_mask,cx,cy,w,Fx,Fy);

        // Bounce-back at solid walls
        bounce_back(Nx,Ny,Q,f,f_new,solid_mask,cx,cy,opp,isperiodic_x,isperiodic_y);

        // Streaming
        streaming(Nx,Ny,Q,f,f_new,solid_mask,cx,cy,isperiodic_x,isperiodic_y);

        // Periodic y BC
        periodic_y(Nx,Ny,Q,0,0,0,f,f_new,rho,u,v,solid_mask,cx,cy,w);

        // Pressure inlet BC
        pressure_inlet(Nx,Ny,Q,0,r0,0,f,f_new,rho,u,v,solid_mask,cx,cy,w);

        // Compute macroscopic quantities
        compute_macroscopic_fields(Nx,Ny,Q,f_new,solid_mask,cx,cy,rho,u,v,Fx,Fy);

        // Compute and save acoustic quantities on the liner surface
        rho_surf = 0.0;
        v_ac = 0.0;
        for (int j = 0; j < Ny; j++) {
            rho_surf += rho[j][i_surf];
            v_ac += u[j][i_surf] * cu;
        }
        v_ac /= Ny;
        rho_surf /= Ny;
        p_ac = c0 * c0 * (rho_surf - rho_0);

        fprintf(fp, "%lf %lf %lf\n", it * dt, p_ac, v_ac);

        // Swap f and f_new
        double (*temp_ptr)[Nx][Q] = f;
        f = f_new;
        f_new = temp_ptr;

        // Write solution (rho, u, v) on a ".vtk" file each save_iter iterations
        if (((it - init_save) % save_iter == 0) & (it >= init_save)) {
            char filename[256];
            snprintf(filename, sizeof(filename),"sol/fields_%gHz_%gPa_%06d.vtk", f_exc, p_a, it);
            write_vtk_binary_2D(filename, Nx, Ny, dx, u, v, rho, cu, crho);
        }

    }

    // Close the tdibc file
    fclose(fp);

    // Free the memory
    free(U_in);
    free(V_in);
    free(rho_in);
    free(solid_mask);
    free(omega_eff);
    free(rho);
    free(u);
    free(v);
    free(f);
    free(f_new);

    return 0;

}