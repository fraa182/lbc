void velocity_outlet(
    int Nx, 
    int Ny,
    int Q, 
    int index,
    double val1, 
    double val2,
    double f[Ny][Nx][Q],
    double f_new[Ny][Nx][Q], 
    double rho[Ny][Nx], 
    double u[Ny][Nx], 
    double v[Ny][Nx], 
    int solid_mask[Ny][Nx], 
    int cx[],
    int cy[],
    double w[]
) {
    double U_out = val1; 
    double V_out = val2;
    int i_out = index;
    
    #pragma omp parallel for
    for (int j = 0; j < Ny; j++) {

        if (solid_mask[j][i_out]) continue;

        double f0 = f_new[j][i_out][0];
        double f2 = f_new[j][i_out][2];
        double f4 = f_new[j][i_out][4];
        double f1 = f_new[j][i_out][1];
        double f5 = f_new[j][i_out][5];
        double f8 = f_new[j][i_out][8];

        double rho_out = (f0 + f2 + f4 + 2.0*(f1 + f5 + f8)) / (1.0 + U_out);

        f_new[j][i_out][3] = f1 - (2.0/3.0)*rho_out*U_out;
        f_new[j][i_out][6] = f8 - 0.5*(f2 - f4) - (1.0/6.0)*rho_out*U_out + (1.0/2.0)*rho_out*V_out;
        f_new[j][i_out][7] = f5 + 0.5*(f2 - f4) - (1.0/6.0)*rho_out*U_out - (1.0/2.0)*rho_out*V_out;

        rho[j][i_out] = rho_out;
        u[j][i_out]   = U_out;
        v[j][i_out]   = V_out;
    }
}
