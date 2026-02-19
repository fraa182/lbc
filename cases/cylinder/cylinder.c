#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include "lbm.h"
#define Q 9
#define save_iter 100

int main(int argc, char *argv[]){

    // Check if the correct number of arguments is provided
    if (argc < 5) {
        fprintf(stderr, "Usage: %s <res> <Re> <tau> <iterations> <use_IBB>\n", argv[0]);
        fprintf(stderr, "Example: %s 10 100 0.9 10000 1\n", argv[0]);
        return 1;
    }

    // UI quantities
    int res = atoi(argv[1]);                      // Resolution (voxels per char length) [-]
    double Reynolds = atof(argv[2]);              // Reynolds number based on U_inf and D = 2R [-]
    double tau = atof(argv[3]);                   // Relaxation time [-]
    int Nt_max = atoi(argv[4]);                   // Maximum number of timesteps [-]
    int use_IBB = atoi(argv[5]);                  // Use IBB (1: yes, 0: no) [-]

    printf("Running simulation with res=%d, Re=%f, tau=%f, for max %d iterations\n", res, Reynolds, tau, Nt_max);

    // Physical parameters (physical units)
    double R = 1.0;                               // Char length [m]
    double Lx = 50*R;                             // Domain x-size [m]
    double Ly = 20*R;                             // Domain y-size [m]
    double xc = Lx/4;                             // Center of the cylinder along x axis [m]
    double yc = Ly/2;                             // Center of the cylinder along y axis [m] 
    double U_inf = 1.0;                           // Inlet axial velocity [m/s]
    double nu = U_inf * (2*R) / Reynolds;         // Kinematic viscosity [m^2/s]
    double dCD_toll = 1e-6;                       // Relative tolerance on CD [-]

    // Lattice domain (lattice units)
    double dx = R / res;                          // Voxel size [m]
    int Nx = Lx / dx;                             // Number of voxels along x
    int Ny = Ly / dx;                             // Number of voxels along y
    double rho_inf = 1.0;                         // Char lattice density

    // Conversion factors (physical to lattice units)
    double crho = rho_inf;                        // Lattice density conversion factor [kg/m^3]
    double cu = 3.0*nu/(dx*(tau - 0.5));          // Lattice velocity conversion factor [m/s]
    double dt = dx / cu;                          // Time step [s]
    double cf = (crho*cu*cu*dx);                  // Lattice force conversion factor [N] -> we are in 2D so tehre is just dx!
    double fref = (0.5*rho_inf*U_inf*U_inf*2*R);  // Reference force to compute CL and CD

    // Forcing (lattice units)
    double Fx = 0.0;
    double Fy = 0.0;

    // Initial conditions (lattice units)
    double U_in = U_inf / cu;                     // Initial horizontal velocity
    double V_in = 0.0;                            // Initial vertical velocity
    double rho_in = rho_inf;                      // Initial density

    // Check lattice Mach number
    double cs = 1.0/sqrt(3.0);
    double Ma = fabs(U_in)/cs;
    if (Ma > 0.2) fprintf(stderr,"Warning: Lattice Ma too high (%.2f), reduce U_inf or increase res or change tau to stay below 0.20\n", Ma);

    // Boundary conditions
    int num_boundaries = 3;
    Boundary boundaries[num_boundaries];

    boundaries[0].apply = periodic_y;

    boundaries[1].index = 0;
    boundaries[1].val1 = U_in;
    boundaries[1].val2 = 0.0;
    boundaries[1].apply = velocity_inlet; 

    boundaries[2].index = Nx - 1;
    boundaries[2].val1 = rho_inf;
    boundaries[2].apply = convective_outlet;

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
    double (*phi)[Nx] = malloc(Ny * sizeof *phi);
    for (int j = 0; j < Ny; j++){
        double y = j*dx;
        for (int i = 0; i < Nx; i++){
            double x = i*dx;
            solid_mask[j][i] = ((x - xc)*(x - xc) + (y - yc)*(y - yc)) <= R*R ? 1 : 0;
            phi[j][i] = sqrt((x - xc)*(x - xc) + (y - yc)*(y - yc)) - R;
        }
    }

    // BGK collision operator coefficient (Omega = omega*(f - f_eq))
    double (*omega_eff)[Nx] = malloc(Ny * sizeof *omega_eff);
    for (int j = 0; j < Ny; j++){
        for (int i = 0; i < Nx; i++){
            omega_eff[j][i] = 1.0/tau;
        }
    }

    // Flow field initialization (in lattice units)
    double (*rho)[Nx] = malloc(Ny * sizeof *rho);
    double (*u)[Nx] = malloc(Ny * sizeof *u);
    double (*v)[Nx] = malloc(Ny * sizeof *v);
    for (int j = 0; j < Ny; j++){
        for (int i = 0; i < Nx; i++){
            u[j][i] = U_in;
            v[j][i] = V_in;
            rho[j][i] = rho_in;
        }
    }
    
    // Particle distribution function initialization
    double (*f)[Nx][Q] = malloc(Ny * sizeof *f);
    double (*f_new)[Nx][Q] = malloc(Ny * sizeof *f_new);
    for (int j = 0; j < Ny; j++){
        for (int i = 0; i < Nx; i++){
            double ux = u[j][i];
            double uy = v[j][i];
            double v_sq = ux*ux + uy*uy;

            for (int k = 0; k < Q; k++){
                double c_dot_u  = cx[k]*ux + cy[k]*uy;
                f[j][i][k] = w[k]*rho[j][i]*(1.0 + 3.0*c_dot_u + 4.5*c_dot_u*c_dot_u - 1.5*v_sq);;
            }
        }
    }
    memcpy(f_new, f, Ny*sizeof(*f));

    // Define and initialize the lift and drag forces on the cylinder surface (physical units)
    double L = 0.0;
    double D = 0.0;
    double CD = 0.0;
    double CD_prec = 1.0;
    double dCD = 1.0;

    // Ensure that the "sol" directory exists and, if not, create it
    ensure_directory_exists("sol");

    char filename_forces[256];
    snprintf(filename_forces, sizeof(filename_forces),"sol/forces_res_%d_Re_%.2f_tau_%.2f.bin", res, Reynolds, tau);

    FILE *fp = fopen(filename_forces, "w");

    struct Record {
        double time;
        double L;
        double D;
    };

    // Temporal loop
    struct Record rec;
    for (int it = 0; it < Nt_max; it++){

        // Compute CD
        CD = D*cf/fref;

        // Print the timestep
        if (it % save_iter == 0) printf("Step %d of %d - CD = %.4f - dCD = %g\n ",it+1, Nt_max, CD, dCD);

        // Execute the streaming and colliding steps
        main_lbm(Nx,Ny,Q,f,f_new,rho,u,v,cx,cy,w,opp,omega_eff,solid_mask,boundaries,num_boundaries,isperiodic_x,isperiodic_y,Fx,Fy,&L,&D,phi,use_IBB);

        // Write forces
        rec.time = it*dt;
        rec.L = L*cf;
        rec.D = D*cf;
        fwrite(&rec, sizeof(struct Record), 1, fp);

        // Swap f and f_new
        double (*temp_ptr)[Nx][Q] = f;
        f = f_new;
        f_new = temp_ptr;

        // Write solution (rho, u, v) on a ".vtk" file each save_iter iterations
        if (it % save_iter == 0){
            char filename[256];
            snprintf(filename, sizeof(filename),"sol/fields_res_%d_Re_%f_tau_%f_%05d.vtk", res, Reynolds, tau, it);
            write_vtk_binary_2D(filename,Nx,Ny,dx,u,v,rho,cu,crho);
        }

        // Check convergence on CD
        dCD = fabs(CD - CD_prec);
        if (dCD <= dCD_toll) {
            printf("Stopping at iteration %d (dCD=%g)\n", it, dCD);
            break;
        } else {
            CD_prec = CD;
        }

    }

    // Close forces file
    fclose(fp);

    // Free the memory
    free(solid_mask);
    free(omega_eff);
    free(rho);
    free(u);
    free(v);
    free(f);
    free(f_new);

    return 0;

}
