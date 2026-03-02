#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include "lbm.h"
#define Q 9
#define save_iter 1

int main(int argc, char *argv[]){

    // Check if the correct number of arguments is provided
    if (argc < 4) {
        fprintf(stderr, "Usage: %s <tau> <# decay times> <initialization iterations>\n", argv[0]);
        return 1;
    }

    // UI quantities
    double tau = atof(argv[1]);                   // Relaxation time [-]
    int Td_sim = atoi(argv[2]);                   // Number of decay times to be simulated [-]
    int init_iter_max = atoi(argv[3]);            // Number of iterations for initialization [-]

    printf("Running simulation with tau = %.2f, Td_sim = %d, init_iter_max = %d\n", tau, Td_sim, init_iter_max);

    // Computational domain (lattice units)
    int Nx = 96;                                  // Number of points along x [-]
    int Ny = 72;                                  // Number of points along y [-]
    double U_inf = 0.03;                          // Initial velocity [m/s]
    double rho_inf = 1.0;                         // Initial density [kg/m^3]
    double nu = (tau - 0.5) / 3;                  // Kinematic viscosity [m^2/s]
    
    // Conversion factors (physical to lattice units)
    double dx = 1.0;                              // Voxel size [m]
    double crho = 1.0;                            // Lattice density conversion factor [kg/m^3]
    double cu = 1.0;                              // Lattice velocity conversion factor [m/s]

    // Calculation of simulation time
    double kx = 2*M_PI / (Nx);                    // Wavenumber along x [1/m]
    double ky = 2*M_PI / (Ny);                    // Wavenumber along y [1/m]
    double decay_time = 1 / (nu*(kx*kx + ky*ky)); // Decay time [s]
    int Nt = (Td_sim*decay_time);                 // Number of time steps [-]

    // Forcing (lattice units)
    double Fx = 0.0;
    double Fy = 0.0;

    // Initial conditions (lattice units)
    double tol_rho = 1e-16;                       // Tolerance for density convergence [-]
    double (*U_in)[Nx] = malloc(Ny * sizeof *U_in);
    double (*V_in)[Nx] = malloc(Ny * sizeof *V_in);
    for (int j = 0; j < Ny; j++){
        for (int i = 0; i < Nx; i++){ 
            U_in[j][i] = -U_inf / cu * sqrt(ky/kx)*cos(kx*i*dx)*sin(ky*j*dx);
            V_in[j][i] = U_inf / cu * sqrt(kx/ky)*sin(kx*i*dx)*cos(ky*j*dx);
        }
    }

    // Boundary conditions
    int use_IBB = 0;                              // Use IBB (1: yes, 0: no) [-]
    int num_boundaries = 2;                       // Number of boundary conditions [-]

    Boundary boundaries[num_boundaries];

    boundaries[0].apply = periodic_x;
    boundaries[1].apply = periodic_y;

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
        for (int i = 0; i < Nx; i++){
            solid_mask[j][i] = 0;
            phi[j][i] = 0.0;
        }
    }

    // BGK collision operator coefficient (Omega = omega*(f - f_eq))
    double (*omega_eff)[Nx] = malloc(Ny * sizeof *omega_eff);
    for (int j = 0; j < Ny; j++){
        for (int i = 0; i < Nx; i++){
            omega_eff[j][i] = 1.0/tau;
        }
    }

    // Flow field and particle distribution function initialization
    double (*rho)[Nx] = malloc(Ny * sizeof *rho);
    double (*u)[Nx] = malloc(Ny * sizeof *u);
    double (*v)[Nx] = malloc(Ny * sizeof *v);
    double (*f)[Nx][Q] = malloc(Ny * sizeof *f);
    double (*f_new)[Nx][Q] = malloc(Ny * sizeof *f_new);
    mei_initialization(Nx, Ny, Q, &f, &f_new, rho_inf, U_in, V_in, rho, u, v, cx, cy, w, omega_eff, isperiodic_x, isperiodic_y, init_iter_max, tol_rho);

    // Define and initialize the lift and drag forces on the surface (physical units)
    double L = 0.0;
    double D = 0.0;

    // Ensure that the "sol" directory exists and, if not, create it
    ensure_directory_exists("sol");

    for (int it = 0; it < Nt; it++){

        // Print the timestep
        if (it % save_iter == 0) printf("Step %d of %d\n ",it+1, Nt);

        // Execute the streaming and colliding steps
        main_lbm(Nx,Ny,Q,f,f_new,rho,u,v,cx,cy,w,opp,omega_eff,solid_mask,boundaries,num_boundaries,isperiodic_x,isperiodic_y,Fx,Fy,&L,&D,phi,use_IBB);

        // Swap f and f_new
        double (*temp_ptr)[Nx][Q] = f;
        f = f_new;
        f_new = temp_ptr;

        // Write solution (rho, u, v) on a ".vtk" file each save_iter iterations
        if (it % save_iter == 0){
            char filename[256];
            snprintf(filename, sizeof(filename),"sol/fields_tau_%f_%05d.vtk", tau, it);
            write_vtk_binary_2D(filename,Nx,Ny,dx,u,v,rho,cu,crho);
        }

    }

    // Free the memory
    free(U_in);
    free(V_in);
    free(solid_mask);
    free(phi);
    free(omega_eff);
    free(rho);
    free(u);
    free(v);
    free(f);
    free(f_new);

    return 0;

}
