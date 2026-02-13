#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
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
    double Reynolds = atof(argv[2]);              // Reynolds number based on U_inf and D = 2R [-]
    int Nt_max = atoi(argv[3]);                   // Maximum number of timesteps [-]

    printf("Running simulation with res=%d, Re=%f for max %d iterations\n", res, Reynolds, Nt_max);

    // Physical parameters
    double R = 0.001;                             // Char length [m]
    double dx = R / res;                          // Voxel size [m]
    double Lx = 0.05;                             // Domain x-size [m]
    double Ly = 0.02;                             // Domain y-size [m]

    // Lattice domain
    int Nx = Lx / dx;                             // Number of voxels along x [-]
    int Ny = Ly / dx;                             // Number of voxels along y [-]
    double U_in = 0.1;                            // Char lattice velocity [m/s]
    double rho0 = 1.0;                            // Char lattice density [kg/m^3]

    // Lattice quantities
    double nu = U_in*(2*R)/Reynolds;              // Kinematic viscosity [m^2/s]
    double crho = rho0;                           // Lattice density conversion factor [kg/m^3]
    double tau = 0.9;                             // Relaxation time in lattice units [-]
    double cu = 3.0*nu/(dx*(tau - 0.5));          // Lattice velocity conversion factor [m/s]

    // Forcing (lattice units)
    double Fx = 0.0;
    double Fy = 0.0;

    // Boundary conditions
    int num_boundaries = 2;
    Boundary boundaries[num_boundaries];

    boundaries[0].index = 0;
    boundaries[0].val1 = U_in;
    boundaries[0].val2 = 0.0;
    boundaries[0].apply = velocity_inlet; 

    boundaries[1].index = Nx - 1;
    boundaries[1].val1 = rho0;
    boundaries[1].apply = convective_outlet;

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

    // Solid mask (staircase approximation)
    int (*solid_mask)[Nx] = malloc(Ny * sizeof *solid_mask);
    //int (*has_wall)[Nx] = malloc(Ny * sizeof *solid_mask);
    //double (*q)[Nx] = malloc(Ny * sizeof *solid_mask);
    for (int j = 0; j < Ny; j++){
        double y = j*dx;
        for (int i = 0; i < Nx; i++){
            double x = i*dx;
            solid_mask[j][i] = ((x - Lx/4)*(x - Lx/4) + (y - Ly/2)*(y - Ly/2)) <= R*R ? 1 : 0;
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
            v[j][i] = 0.0;
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
                f[j][i][k] = w[k]*rho[j][i]*(1.0 + 3.0*cu + 4.5*cu*cu - 1.5*v_sq);;
            }
        }
    }
    memcpy(f_new, f, Ny*sizeof(*f));

    // Temporal loop
    for (int it = 0; it < Nt_max; it++){
        // Print the timestep
        if (it % save_iter == 0) printf("Step %d of %d\n",it+1,Nt_max);

        // Execute the streaming and colliding steps
        main_lbm(Nx,Ny,Q,f,f_new,rho,u,v,cx,cy,w,opp,omega_eff,solid_mask,boundaries,num_boundaries,isperiodic_x,isperiodic_y,Fx,Fy);

        // Swap f and f_new
        double (*temp_ptr)[Nx][Q] = f;
        f = f_new;
        f_new = temp_ptr;

        // Ensure that the "sol" directory exists and, if not, create it
        ensure_directory_exists("sol");

        // Write solution (rho, u, v) on a ".vtk" file each save_iter iterations
        if (it % save_iter == 0){
            char filename[256];
            snprintf(filename, sizeof(filename),"sol/fields_%05d.vtk", it);
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
