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
    if (argc < 4) {
        fprintf(stderr, "Usage: %s <res> <tau> <max iterations>\n", argv[0]);
        fprintf(stderr, "Example: %s 10 0.90 100000\n", argv[0]);
        return 1;
    }

    // UI quantities
    int res = atoi(argv[1]);                      // Resolution (voxels per char length) [-]
    double tau = atof(argv[2]);                   // Relaxation time in lattice units [-]
    int Nt_max = atoi(argv[3]);                   // Maximum number of timesteps [-]

    printf("Running simulation with res=%d, tau=%f for max %d iterations\n", res, tau, Nt_max);

    // Safety check for tau (LBM stability)
    if (tau <= 0.5) {
        fprintf(stderr, "Error: tau must be greater than 0.5 for stability.\n");
        return 1;
    }

    // Relative tolerance on centreline velocity convergence
    double dvmax_toll = 1e-16;

    // Lattice domain (we work in lattice units here so dx = 1 and cu = 1, they are just for export)
    int Nx = 50;                                  // Number of voxels along x [-]
    int Ny = res;                                 // Number of voxels along y [-]
    double dx = 1.0;                              // Voxel size [m]
    double rho0 = 1.0;                            // Char lattice density [kg/m^3]
    double crho = rho0;                           // Lattice density conversion factor [kg/m^3]
    double cu = 1.0;                              // Lattice velocity conversion factor [m/s]

    // Forcing (lattice units)
    double Fx = 0.0;
    double Fy = 0.0;

    // Boundary conditions
    int num_boundaries = 2;                       // Number of boundaries [-]
    Boundary boundaries[num_boundaries];          // Define the Boundary struct

    boundaries[0].apply = periodic_x;             // BC type (1)
    
    boundaries[1].apply = velocity_top;           // BC type (2)
    boundaries[1].index = Ny - 1;                 // BC index (2)
    boundaries[1].val1 = 0.1;                     // BC value for U (2)
    boundaries[1].val2 = 0.0;                     // BC value for V (2)
    
    // Check if there is a periodic BC along x and/or y
    int isperiodic_x = 0;
    int isperiodic_y = 0;
    for (int i = 0; i < num_boundaries; i++) {
        if (boundaries[i].apply == periodic_x) isperiodic_x = 1;
        if (boundaries[i].apply == periodic_y) isperiodic_y = 1;
    }

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
            solid_mask[j][i] = (j == 0) ? 1 : 0;
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
            u[j][i] = 0.0;
            v[j][i] = 0.0;
            rho[j][i] = rho0;
        }
    }
    
    // Particle distribution function initialization
    double (*f)[Nx][Q] = malloc(Ny * sizeof *f);
    double (*f_new)[Nx][Q] = malloc(Ny * sizeof *f_new);
    for (int j = 0; j < Ny; j++){
        for (int i = 0; i < Nx; i++){
            for (int k = 0; k < Q; k++){
                f[j][i][k] = w[k]*rho0;
            }
        }
    }
    memcpy(f_new, f, Ny*sizeof(*f_new));

    // Initialize velocity difference and centreline axial velocity at previous step
    double dv = 1;
    double ux_center_prec = 1;

    // Temporal loop
    for (int it = 0; it < Nt_max; it++){
        // Execute the streaming and colliding steps
        main_lbm(Nx,Ny,Q,f,f_new,rho,u,v,cx,cy,w,opp,omega_eff,solid_mask,boundaries,num_boundaries,isperiodic_x,isperiodic_y,Fx,Fy);

        // Swap f and f_new
        double (*temp_ptr)[Nx][Q] = f;
        f = f_new;
        f_new = temp_ptr;

        // Check the maximum lattice velocity
        double ux_ref = 0.0;
        for (int i = 0; i < Nx; i++){
            ux_ref += u[Ny-2][i];
        }
        ux_ref /= (double)Nx;

        double umax = 0.0;
        for (int j = 0; j < Ny; j++){
            for (int i = 0; i < Nx; i++){
                double au = fabs(u[j][i]);
                if (au > umax) umax = au;
            }
        }

        if (umax > 0.2) { 
            fprintf(stderr, "Error: max |u| too high: %g\n", umax);
            return 1;
        }

        // Print the timestep
        if (it % save_iter == 0) printf("Step %d of %d - ref velocity: %f\n", it+1, Nt_max, ux_ref);

        // Ensure that the "sol" directory exists and, if not, create it
        ensure_directory_exists("sol");

        // Write solution (rho, u, v) on a ".vtk" file each save_iter iterations
        if (it % save_iter == 0){
            char filename[256];
            snprintf(filename, sizeof(filename),"sol/fields_res_%d_tau_%.2f_%05d.vtk", res, tau, it);
            write_vtk_binary_2D(filename,Nx,Ny,dx,u,v,rho,cu,crho);
        }

        // Exit from the loop if the maximum velocity does not change from the previous iteration (within tolerance)
        double denom = fabs(ux_ref);
        if (denom > 1e-12){
            dv = fabs(ux_ref - ux_center_prec)/denom;
        } else {
            dv = fabs(ux_ref - ux_center_prec);
        }
        
        ux_center_prec = ux_ref;

        if (dv <= dvmax_toll) {
            printf("Stopping at iteration %d (dv=%g)\n", it, dv);
            
            char filename[256];
            snprintf(filename, sizeof(filename),"sol/fields_res_%d_tau_%.2f_%05d.vtk", res, tau, it);
            write_vtk_binary_2D(filename,Nx,Ny,dx,u,v,rho,cu,crho);

            break;
        }
    }

    // Let everyone know that the tolerance criterion was not met
    if (dv > dvmax_toll) printf("Tolerance criterion not reached. Increase max iterations.\n");

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