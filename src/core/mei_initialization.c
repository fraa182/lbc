#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

void mei_initialization(
    int Nx,
    int Ny, 
    int Q,
    double (**pf)[Nx][Q], 
    double (**pf_new)[Nx][Q], 
    double rho_inf,
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
    double tol_rho
){

    printf("Initialization of the populations... ");

    // Initialize velocity and density fields
    for (int j = 0; j < Ny; j++){
        for (int i = 0; i < Nx; i++){
            u[j][i] = U_in[j][i];
            v[j][i] = V_in[j][i];
            rho[j][i] = rho_inf;
        }
    }

    // Define local pointers of f and f_new
    double (*f)[Nx][Q] = *pf;
    double (*f_new)[Nx][Q] = *pf_new;

    // Initialize particle distribution function (incompressible equilibrium)
    for (int j = 0; j < Ny; j++){
        for (int i = 0; i < Nx; i++){
            double ux = u[j][i];
            double uy = v[j][i];
            double v_sq = ux*ux + uy*uy;

            for (int k = 0; k < Q; k++){
                double cu  = cx[k]*ux + cy[k]*uy;
                f[j][i][k] = w[k]*rho[j][i] + w[k]*rho_inf*(3.0*cu + 4.5*cu*cu - 1.5*v_sq);;
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
                    double feq = w[k]*rho[j][i] + w[k]*rho_inf*(3.0*cu + 4.5*cu*cu - 1.5*v_sq);

                    f[j][i][k] += -omega_eff[j][i] * (f[j][i][k] - feq);
                }
            }
        }

        // Perform streaming
        #pragma omp parallel for
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {

                for (int k = 0; k < Q; k++) {
                    int i_src = i - cx[k];
                    int j_src = j - cy[k];

                    // If we're using the periodic BC along x, wrap x
                    if (isperiodic_x){
                        if (i_src < 0) i_src += Nx;
                        else if (i_src >= Nx) i_src -= Nx;
                    }

                    // If we're using the periodic BC along y, wrap y
                    if (isperiodic_y){
                        if (j_src < 0) j_src += Ny;
                        else if (j_src >= Ny) j_src -= Ny;
                    }

                    // If non-periodic and out of bounds: leave it for BCs to set later
                    if (i_src < 0 || i_src >= Nx || j_src < 0 || j_src >= Ny) continue;

                    f_new[j][i][k] = f[j_src][i_src][k];
                }
            }
        }

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

    // Write back the final pointers to caller
    *pf = f;
    *pf_new = f_new;

}