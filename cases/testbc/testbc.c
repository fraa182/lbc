#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include "lbm.h"
#define Q 9
#define save_iter 10
#define disp_iter 100

// Main function
int main(int argc, char *argv[]){

    // Computational domain (physical units)
    int Nt = 500;                                 // Number of time steps [-]
    double Lx = 0.010;                            // Domain length along x [m]
    double Ly = 0.010;                            // Domain length along y [m]
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

    // Forcing (lattice units)
    double Fx = 0.0;                              // Body force x component [m / s^2]
    double Fy = 0.0;                              // Body force y component [m / s^2]

    // Initial conditions (lattice units)
    int init_iter_max = 100;                      // Maximum number of iterations in Mei's algorithm [-]
    double tol_rho = 1e-10;                       // Tolerance for density convergence in Mei's algorithm [-]

    // Compute derivate quantities
    int Nx = ceil(Lx / dx);                       // Number of points along x [-]
    int Ny = ceil(Ly / dx);                       // Number of points along y [-]

    // Boundary conditions
    int use_IBB = 0;                              // Use IBB (1: yes, 0: no) [-]
    int num_boundaries = 4;                       // Number of boundary conditions [-]

    Boundary boundaries[num_boundaries];

    double dmy_inout[Ny];
    for (int j = 0; j < Ny; j++) {
        dmy_inout[j] = 0.0;
    } 

    double dmy_topbot[Nx];
    for (int i = 0; i < Nx; i++) {
        dmy_topbot[i] = 0.0;
    }

    boundaries[0].apply = velocity_inlet_regularized;
    boundaries[0].index = 0;
    boundaries[0].val1 = dmy_inout;
    boundaries[0].val2 = dmy_inout;

    boundaries[1].apply = velocity_outlet_regularized;
    boundaries[1].index = Nx-1;
    boundaries[1].val1 = dmy_inout;
    boundaries[1].val2 = dmy_inout;

    boundaries[2].apply = velocity_bottom_regularized;
    boundaries[2].index = 0;
    boundaries[2].val1 = dmy_topbot;
    boundaries[2].val2 = dmy_topbot;

    boundaries[3].apply = velocity_top_regularized;
    boundaries[3].index = Ny-1;
    boundaries[3].val1 = dmy_topbot;
    boundaries[3].val2 = dmy_topbot;

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
            if (i <= Nx / 4) {
                omega_eff[j][i] = 1.0 / 0.9;
            } else {
                omega_eff[j][i] = 1.0 / tau;
            }
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
            rho_in[j][i] = rho_0;

            U_in[j][i] = 0.0;
            V_in[j][i] = 0.0;
        }
    }
    mei_initialization(Nx, Ny, Q, f, f_new, rho_in, U_in, V_in, rho, u, v, cx, cy, w, omega_eff, isperiodic_x, isperiodic_y, init_iter_max, tol_rho, solid_mask);

    // Define and initialize the lift and drag forces on the surface (physical units)
    double L = 0.0;
    double D = 0.0;

    double v_bc = 0.0;
    double vBC_inout[Ny];
    double vBC_topbot[Nx];

    // Ensure that the "sol" directory exists and, if not, create it
    ensure_directory_exists("sol");

    // Main LBM loop
    for (int it = 0; it < Nt; it++){

        // Print the timestep
        if (it % disp_iter == 0) printf("Step %d of %d\n ", it, Nt);

        // Define boundary velocity (physical units)
        v_bc = 2 * sin(2 * M_PI * it / 50.0);

        // Fill the local values of velocity (lattice units)
        for (int j = 0; j < Ny; j++){
            vBC_inout[j] = (v_bc / cu) * sin(2 * M_PI * 2 * j * dx / Ly);
        }

        for (int i = 0; i < Nx; i++){
            vBC_topbot[i] = (v_bc / cu) * sin(2 * M_PI * 2 * i * dx / Lx);
        }

        // Apply velocity BC (lattice units)
        boundaries[0].val1 = vBC_inout;
        boundaries[1].val1 = vBC_inout;
        boundaries[2].val2 = vBC_topbot;
        boundaries[3].val2 = vBC_topbot;

        // Execute the streaming and colliding steps
        main_lbm(Nx,Ny,Q,f,f_new,rho,u,v,cx,cy,w,opp,omega_eff,solid_mask,boundaries,num_boundaries,isperiodic_x,isperiodic_y,Fx,Fy,&L,&D,phi_solid,use_IBB);

        // Swap f and f_new
        double (*temp_ptr)[Nx][Q] = f;
        f = f_new;
        f_new = temp_ptr;

        // Write solution (rho, u, v) on a ".vtk" file each save_iter iterations
        if (it % save_iter == 0){
            char filename[256];
            snprintf(filename, sizeof(filename),"sol/fields_%05d.vtk", it);
            write_vtk_binary_2D(filename, Nx, Ny, dx, u, v, rho, cu, crho);
        }

    }

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
