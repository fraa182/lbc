void velocity_inlet(
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
    double U_in = val1;
    double V_in = val2;
    int i_in = index;
    
    #pragma omp parallel for
    for (int j = 0; j < Ny; j++) {

        if (solid_mask[j][i_in]) continue;

        double f0 = f_new[j][i_in][0];
        double f2 = f_new[j][i_in][2];
        double f4 = f_new[j][i_in][4];
        double f3 = f_new[j][i_in][3];
        double f6 = f_new[j][i_in][6];
        double f7 = f_new[j][i_in][7];

        double rho_in = (f0 + f2 + f4 + 2.0*(f3 + f6 + f7)) / (1.0 - U_in);

        f_new[j][i_in][1] = f3 + (2.0/3.0)*rho_in*U_in;
        f_new[j][i_in][5] = f7 + (1.0/6.0)*rho_in*U_in - 0.5*(f2 - f4);
        f_new[j][i_in][8] = f6 + (1.0/6.0)*rho_in*U_in + 0.5*(f2 - f4);

        rho[j][i_in] = rho_in;
        u[j][i_in]   = U_in;
        v[j][i_in]   = V_in;
    }
}
