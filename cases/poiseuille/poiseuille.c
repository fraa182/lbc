#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../../src/include/lbm.h"
#define Q 9
#define save_iter 100

int main(int argc, char *argv[]){

    // Check if the correct number of arguments is provided
    if (argc < 4) {
        fprintf(stderr, "Usage: %s <res> <tau> <simulation time>\n", argv[0]);
        fprintf(stderr, "Example: %s 10 0.90 25\n", argv[0]);
        return 1;
    }

    // UI quantities
    int res = atoi(argv[1]);                      // Resolution (voxels per char length) [-]
    double tau = atof(argv[2]);                   // Relaxation time in lattice units [-]
    double tf = atof(argv[3]);                    // Simulation time [s]

    // Safety check for tau (LBM stability)
    if (tau <= 0.5) {
        fprintf(stderr, "Error: tau must be greater than 0.5 for stability.\n");
        return 1;
    }

    printf("Running simulation with res=%d, tau=%.3f for %.2f seconds\n", res, tau, tf);

    // Physical parameters
    
    double nu = 1.5e-5;                           // Kinematic viscosity [m^2/s]
    double Lx = 0.100;                            // Domain x-size [m]
    double Ly = 0.025;                            // Domain y-size (char length) [m]

    // Lattice domain
    double dx = Ly / res;                         // Voxel size [m]
    int Nx = Lx / dx;                             // Number of voxels along x [-]
    int Ny = Ly / dx;                             // Number of voxels along y [-]
    double rho0 = 1.0;                            // Char lattice density [kg/m^3]

    // Forcing
    double Fx = 2.0e-5;
    double Fy = 0.0;

    // Boundary conditions
    int num_boundaries = 1;                       // Number of boundaries [-]
    Boundary boundaries[num_boundaries];          // Define the Boundary struct
    boundaries[0].apply = periodic_x;             // BC type

    // Lattice quantities
    double ct = (tau - 0.5)*(dx*dx)/(3*nu);       // Timestep [s]
    double crho = rho0;                           // Lattice density conversion factor [kg/m^3]
    double cu = dx / ct;                          // Lattice velocity conversion factor [m/s]
    int Nt = tf / ct;                             // Number of timesteps [-]

    // Lattice velocity x and y components (D2Q9)
    int cx[Q] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
    int cy[Q] = {0, 0, 1, 0, -1, 1, 1, -1, -1};

    // Lattice weights (D2Q9)
    double w[Q] = {4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};

    // Indices of opposite directions for D2Q9 lattice
    int opp[Q] = {0, 3, 4, 1, 2, 7, 8, 5, 6};

    // Solid mask
    int (*solid_mask)[Nx] = malloc(Ny * sizeof *solid_mask);
    for (int j = 0; j < Ny; j++){
        for (int i = 0; i < Nx; i++){
            solid_mask[j][i] = ((j == 0) | (j == Ny-1)) ? 1 : 0;
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
            u[j][i] = 0;
            v[j][i] = 0;
            rho[j][i] = 1.0;
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
                double cu  = cx[k]*ux + cy[k]*uy;
                f[j][i][k] = w[k]*rho0*(1.0 + 3.0*cu + 4.5*cu*cu - 1.5*v_sq);;
            }
        }
    }
    memcpy(f_new, f, Ny*sizeof(*f));

    // Temporal loop
    for (int it = 0; it < Nt; it++){
        // Print the timestep
        if (it % save_iter == 0) printf("Step %d of %d\n",it+1,Nt);

        // Execute the streaming and colliding steps
        main_lbm(Nx,Ny,Q,f,f_new,rho,u,v,cx,cy,w,opp,omega_eff,solid_mask,boundaries,num_boundaries,Fx,Fy);

        // Swap f and f_new
        double (*temp_ptr)[Nx][Q] = f;
        f = f_new;
        f_new = temp_ptr;

        // Ensure that the "sol" directory exists and, if not, create it
        ensure_directory_exists("sol");

        // Write solution (rho, u, v) on a ".vtk" file each save_iter iterations
        if (it % save_iter == 0){
            char filename[256];
            snprintf(filename, sizeof(filename),"sol/fields_res_%d_%05d.vtk", res, it);
            write_vtk_binary_2D(filename,Nx,Ny,dx,u,v,rho,cu,crho);
        }
    }

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