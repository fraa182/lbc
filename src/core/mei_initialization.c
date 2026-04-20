#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "lbm.h"

void mei_initialization(
    int Nx,
    int Ny, 
    int Q,
    double f[Ny][Nx][Q], 
    double f_new[Ny][Nx][Q], 
    double rho_in[Ny][Nx],
    double U_in[Ny][Nx], 
    double V_in[Ny][Nx],
    double rho[Ny][Nx], 
    double u[Ny][Nx], 
    double v[Ny][Nx], 
    int cx[], 
    int cy[], 
    double w[], 
    double omega_eff[Ny][Nx], 
    int isperiodic_x,
    int isperiodic_y,
    int init_iter_max,
    double tol_rho,
    int solid_mask[Ny][Nx]
){

    printf("Initialization of the populations... ");

    // Initialize velocity and density fields
    for (int j = 0; j < Ny; j++){
        for (int i = 0; i < Nx; i++){
            u[j][i] = U_in[j][i];
            v[j][i] = V_in[j][i];
            rho[j][i] = rho_in[j][i];
        }
    }

    // Compute mean density for incompressible equilibrium population
    double rho_0 = 0.0;
    for (int j = 0; j < Ny; j++){
        for (int i = 0; i < Nx; i++){
            rho_0 += rho[j][i];
        }
    }
    rho_0 /= (Nx * Ny);

    // Initialize particle distribution function (incompressible equilibrium)
    #pragma omp parallel for
    for (int j = 0; j < Ny; j++){
        for (int i = 0; i < Nx; i++){
            double ux = u[j][i];
            double uy = v[j][i];
            double v_sq = ux*ux + uy*uy;

            for (int k = 0; k < Q; k++){
                double cu  = cx[k]*ux + cy[k]*uy;
                f[j][i][k] = w[k]*rho[j][i] + w[k]*rho_0*(3.0*cu + 4.5*cu*cu - 1.5*v_sq);
            }
        }
    }

    // Temporary array for previous rho to test convergence
    double (*rho_prev)[Nx] = malloc(Ny * sizeof *rho_prev);
    double max_drho = 0.0;

    // Iterate until density convergence or max iterations
    for (int init_iter = 0; init_iter < init_iter_max; init_iter++) {

        // Save previous rho
        memcpy(rho_prev, rho, Ny * sizeof(rho[0]));

        // Compute density
        #pragma omp parallel for
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {

                double s  = 0.0;

                for (int k = 0; k < Q; k++) {
                    double fk = f[j][i][k];
                    s  += fk;
                }

                rho[j][i] = s;
            }
        }

        // Perform incompressible collision
        #pragma omp parallel for
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {

                double ux = u[j][i];
                double uy = v[j][i];

                double v_sq = ux*ux + uy*uy;

                for (int k = 0; k < Q; k++) {
                    double cu  = cx[k]*ux + cy[k]*uy;
                    double feq = w[k]*rho[j][i] + w[k]*rho_0*(3.0*cu + 4.5*cu*cu - 1.5*v_sq);

                    f[j][i][k] += -omega_eff[j][i] * (f[j][i][k] - feq);
                }
            }
        }

        // Perform streaming
        streaming(Nx,Ny,Q,f,f_new,solid_mask,cx,cy,isperiodic_x,isperiodic_y);

        // Swap f and f_new
        double (*temp_ptr)[Nx][Q] = f;
        f = f_new;
        f_new = temp_ptr;

        // Check if the density change from the previous iteration is below the tolerance
        max_drho = 0.0;

        #pragma omp parallel for reduction(max:max_drho)
        for (int j = 0; j < Ny; j++){
            for (int i = 0; i < Nx; i++){
                const double dr = fabs(rho[j][i] - rho_prev[j][i]);
                if (dr > max_drho) max_drho = dr;
            }
        }

        if (max_drho <= tol_rho){
            printf("Initialization converged at iter %d with max |drho| = %g\n", init_iter, max_drho);
            break;
        }

    }

    // If the initialization did not converge, let the user know
    if (max_drho > tol_rho) {
        printf("Initialization not converged after %d iterations (max |drho| = %g)\n", init_iter_max, max_drho);
    }

    // Free the memory
    free(rho_prev);
}