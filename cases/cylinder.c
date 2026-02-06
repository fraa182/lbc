#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../src/include/lbm.h"
#define Q 9
#define save_iter 100

int main(){

    // Physical parameters
    int res = 10;                                // Resolution (voxels per char length) [-]
    double R = 0.001;                           // Char length [m]
    double vx_size = R/res;                     // Voxel size [m]

    double U_phys = 10;                         // Physical velocity [m/s]
    double nu_phys = 0.000015;                  // Kinematic viscosity [m^2/s]
    double Lx = 0.05;                           // Domain x-size [m]
    double Ly = 0.02;                           // Domain y-size [m]
    double t_final = 0.0025;                     // Simulation physical time [s]

    // Lattice domain
    int Nx = Lx / vx_size;                     // Number of voxels along x [-]
    int Ny = Ly / vx_size;                     // Number of voxels along y [-]

    // Boundary conditions (in lattice units)
    double U_in = 0.1;                          // Inlet velocity [-]
    double rho0 = 1.0;                          // Outlet density [-]
    int i_in = 0;                               // Index of inlet [-]
    int i_out = Nx - 1;                         // Index of outlet [-]

    // Lattice quantities
    double cu = U_phys/U_in;                    // Lattice velocity conversion factor [m/s]
    double ct = vx_size/cu;                     // Timestep [s]
    int Nt = t_final/ct;                        // Number of timesteps [-]

    double nu_lattice = nu_phys/(cu*vx_size);   // Lattice viscosity [-]
    double tau_base = 0.5 + 3*nu_lattice;       // Relaxation time in lattice units [-]

    // Lattice velocity x and y components
    int cx[Q] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
    int cy[Q] = {0, 0, 1, 0, -1, 1, 1, -1, -1};

    // Lattice weights (D2Q9)
    double w[Q] = {4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};

    // Indices of opposite directions for [0,...,8]
    int opp[Q] = {0, 3, 4, 1, 2, 7, 8, 5, 6};

    // Solid mask
    int (*solid_mask)[Nx] = malloc(Ny * sizeof *solid_mask);
    for (int j = 0; j < Ny; j++){
        double y = j*vx_size;
        for (int i = 0; i < Nx; i++){
            double x = i*vx_size;
            solid_mask[j][i] = ((x - Lx/2)*(x - Lx/2) + (y - Ly/2)*(y - Ly/2)) <= R*R ? 1 : 0;
        }
    }

    // BGK collision operator coefficient (Omega = omega*(f - f_eq))
    double (*omega_eff)[Nx] = malloc(Ny * sizeof *omega_eff);
    for (int j = 0; j < Ny; j++){
        for (int i = 0; i < Nx; i++){
            omega_eff[j][i] = 1.0/tau_base;
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
            for (int k = 0; k < Q; k++){
                f[j][i][k] = w[k]*rho0;
                f_new[j][i][k] = f[j][i][k];
            }
        }
    }

    // Temporal loop
    for (int it = 0; it < Nt; it++){
        // Print the timestep
        printf("Step %d of %d\n",it+1,Nt);

        // Execute the streaming and colliding steps
        main_lbm(Nx,Ny,f,f_new,rho,u,v,cx,cy,w,opp,omega_eff,solid_mask,i_in,U_in,i_out,rho0);

        // Swap f and f_new
        for (int j = 0; j < Ny; j++){
            for (int i = 0; i < Nx; i++){
                for (int k = 0; k < Q; k++) {
                    double tmp = f[j][i][k];
                    f[j][i][k] = f_new[j][i][k];
                    f_new[j][i][k] = tmp;
                }
            }
        }

        // Ensure that the "sol" directory exists and, if not, create it
        ensure_directory_exists("sol");

        // Write solution (rho, u, v) on a ".vtk" file each save_iter iterations
        if (it % save_iter == 0){
            char filename[256];
            snprintf(filename, sizeof(filename),"sol/fields_%05d.vtk", it);
            write_vtk_binary_2D(filename,Nx,Ny,vx_size,u,v,rho);
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
