#include <stdio.h>
#define Q 9
#define myzero 1e-16

void lbm_step(int Nx, int Ny, double f[Ny][Nx][Q], double f_new[Ny][Nx][Q], double rho[Ny][Nx], double u[Ny][Nx], double v[Ny][Nx], int cx[], int cy[], double w[], int opp[], double omega_eff[Ny][Nx], int solid_mask[Ny][Nx], int i_in, double U_in, int i_out, double rho_out){

    // Collision step
    for (int i = 0; i < Nx; i++){
        for (int j = 0; j < Ny; j++){
            if (!solid_mask[j][i]){
                double v_sq = u[j][i]*u[j][i] + v[j][i]*v[j][i];
                for (int k = 0; k < Q; k++){
                    double cu = cx[k]*u[j][i] + cy[k]*v[j][i];
                    double feq = w[k]*rho[j][i]*(1.0 + 3.0*cu + 4.5*cu*cu - 1.5*v_sq);
                    f[j][i][k] += -omega_eff[j][i] * (f[j][i][k] - feq);
                }
            }
        }
    }

    // Streaming step
    for (int k = 0; k < Q; k++){
        for (int i = 0; i < Nx; i++){
            int i_src = i - cx[k];
            if (i_src < 0) i_src += Nx;
            if (i_src >= Nx) i_src -= Nx;
            for (int j = 0; j < Ny; j++){
                int j_src = j - cy[k];
                if (j_src < 0) j_src += Ny;
                if (j_src >= Ny) j_src -= Ny;
                f_new[j][i][k] = f[j_src][i_src][k];
            }
        }
    }

    // Bounce-back boundary condition at solid walls
    for (int i = 0; i < Nx; i++){
        for (int j = 0; j < Ny; j++){
            if (solid_mask[j][i]){
                double tmp[Q];
                for (int k = 0; k < Q; k++) tmp[k] = f_new[j][i][k];
                for (int k = 0; k < Q; k++) f_new[j][i][k] = tmp[opp[k]];
            }
        }
    }

    // Velocity inlet
    for (int j = 0; j < Ny; j++){
        if (!solid_mask[j][i_in]){
            u[j][i_in]  = U_in;
            v[j][i_in]  = 0.0;
            double v_sq_in = u[j][i_in]*u[j][i_in] + v[j][i_in]*v[j][i_in];
            for (int k = 0; k < Q; k++){
                double cu = cx[k]*u[j][i_in] + cy[k]*v[j][i_in];
                f_new[j][i_in][k] = w[k]*rho_out*(1.0 + 3.0*cu + 4.5*cu*cu - 1.5*v_sq_in);
            }
        }
    }

    // Pressure outlet
    int i_int = i_out > 0 ? i_out - 1 : i_out + 1;
    for (int j = 0; j < Ny; j++){
        if (!solid_mask[j][i_out]){
            double s = 0.0;
            double sx = 0.0;
            double sy = 0.0;
            for (int k = 0; k < Q; k++){
                s  += f_new[j][i_int][k];
                sx += f_new[j][i_int][k] * cx[k];
                sy += f_new[j][i_int][k] * cy[k];
            }
            if (s < myzero){
                s = myzero;
            }
            double u_int = sx / s;
            double v_int = sy / s;

            double v_sq_out = u_int*u_int + v_int*v_int;
            rho[j][i_out] = rho_out;
            u[j][i_out]  = u_int;
            v[j][i_out]  = v_int;
            for (int k = 0; k < Q; k++){
                double cu = cx[k]*u_int + cy[k]*v_int;
                f_new[j][i_out][k] = w[k]*rho_out*(1.0 + 3.0*cu + 4.5*cu*cu - 1.5*v_sq_out);
            }
        }
    }
    
    // Compute macroscopic quantities
    for (int i = 0; i < Nx; i++){
        for (int j = 0; j < Ny; j++){
            double s = 0.0;
            double sx = 0.0;
            double sy = 0.0;
            for (int k = 0; k < Q; k++){
                s  += f_new[j][i][k];
                sx += f_new[j][i][k] * cx[k];
                sy += f_new[j][i][k] * cy[k];
            }
            if (s < myzero){
                s = myzero;
            }
            rho[j][i] = s;
            if (!solid_mask[j][i]){
                u[j][i]  = sx / s;
                v[j][i]  = sy / s;
            }
            else{
                u[j][i] = 0.0;
                v[j][i] = 0.0;
            }
        }
    }

}
