#ifndef MAIN_H
#define MAIN_H
#define Q 9

void main_lbm(
    int Nx,
    int Ny, 
    double f[Ny][Nx][Q], 
    double f_new[Ny][Nx][Q], 
    double rho[Ny][Nx], 
    double u[Ny][Nx], 
    double v[Ny][Nx], 
    int cx[], 
    int cy[], 
    double w[], 
    int opp[], 
    double omega_eff[Ny][Nx], 
    int solid_mask[Ny][Nx], 
    int i_in, 
    double U_in, 
    int i_out, 
    double rho_out
);

#endif
