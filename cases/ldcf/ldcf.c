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
        fprintf(stderr, "Usage: %s <res> <Re> <tau> <max_time>\n", argv[0]);
        fprintf(stderr, "Example: %s 20 100 0.9 10\n", argv[0]);
        return 1;
    }

    // UI quantities
    int res = atoi(argv[1]);                      // Resolution (voxels per char length) [-]
    double Reynolds = atof(argv[2]);              // Reynolds number based on U_inf and D = 2R [-]
    double tau = atof(argv[3]);                   // Relaxation time [-]
    double T_max = atoi(argv[4]);                 // Maximum time [s]

    // Physical parameters (physical units)
    double L = 1.0;                               // Char length [m]
    double U_inf = 1.0;                           // Inlet axial velocity [m/s]
    double nu = U_inf * L / Reynolds;             // Kinematic viscosity [m^2/s]

    // Lattice domain (lattice units)
    double dx = L / res ;                         // Voxel size [m]                       
    int Nx = L / dx;                              // Number of voxels along x
    int Ny = L / dx;                              // Number of voxels along y
    double rho_inf = 1.0;                         // Char lattice density

    // Conversion factors (physical to lattice units)
    double crho = rho_inf;                        // Lattice density conversion factor [kg/m^3]
    double cu = 3.0*nu/(dx*(tau - 0.5));          // Lattice velocity conversion factor [m/s]

    // Max iterations
    int Nt_max = T_max / (dx/cu);

    // Forcing (lattice units)
    double Fx = 0.0;
    double Fy = 0.0;

    // Initial conditions (lattice units)
    double U_in = 0.0;                            // Initial horizontal velocity
    double V_in = 0.0;                            // Initial vertical velocity
    double rho_in = rho_inf;                      // Initial density

    printf("Running simulation with res=%d, Re=%f, tau=%f, for max %d iterations (dt = %f s)\n", res, Reynolds, tau, Nt_max, dx/cu);

    // Check lattice Mach number
    double cs = 1.0/sqrt(3.0);
    double Ma = fabs(U_inf)/(cu*cs);
    if (Ma > 0.2) fprintf(stderr,"Warning: Lattice Ma too high (%.2f), reduce U_inf or increase res or change tau to stay below 0.20\n", Ma);

    // Boundary conditions
    int num_boundaries = 1;
    Boundary boundaries[num_boundaries];

    boundaries[0].index = Ny - 1;
    boundaries[0].val1 = U_inf / cu;
    boundaries[0].val2 = 0.0;
    boundaries[0].apply = velocity_top; 

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

    // Solid mask
    int (*solid_mask)[Nx] = malloc(Ny * sizeof *solid_mask);
    for (int j = 0; j < Ny; j++){
        for (int i = 0; i < Nx; i++){
            solid_mask[j][i] = (j == 0) | (i == 0) | (i == Nx - 1) ? 1 : 0;
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

    int use_IBB = 0;
    double FL = 0.0;
    double FD = 0.0;
    double (*phi)[Nx] = malloc(Ny * sizeof *phi);

    // Ensure that the "sol" directory exists and, if not, create it
    ensure_directory_exists("sol");

    for (int it = 0; it < Nt_max; it++){

        // Print the timestep
        if (it % save_iter == 0) printf("Step %d of %d\n ",it+1, Nt_max);

        // Execute the streaming and colliding steps
        main_lbm(Nx,Ny,Q,f,f_new,rho,u,v,cx,cy,w,opp,omega_eff,solid_mask,boundaries,num_boundaries,isperiodic_x,isperiodic_y,Fx,Fy,&FL,&FD,phi,use_IBB);

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
