void velocity_bottom(
    int Nx, 
    int Ny,
    int Q, 
    int index,
    const double *val1, 
    const double *val2,
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
    int j_bot = index;
    
    #pragma omp parallel for
    for (int i = 0; i < Nx; i++) {

        if (solid_mask[j_bot][i]) continue;

        double U_bot = val1[i];
        double V_bot = val2[i];

        double f0 = f_new[j_bot][i][0];
        double f1 = f_new[j_bot][i][1];
        double f3 = f_new[j_bot][i][3]; 
        double f4 = f_new[j_bot][i][4];
        double f7 = f_new[j_bot][i][7];
        double f8 = f_new[j_bot][i][8];

        double rho_bot = (f0 + f1 + f3 + 2*(f4 + f7 + f8))/(1 - V_bot);

        f_new[j_bot][i][2] = f4 + 2.0/3.0*rho_bot*V_bot;
        f_new[j_bot][i][5] = f7 - (f1 - f3)/2.0 + 0.5*rho_bot*U_bot + 1.0/6.0*rho_bot*V_bot;
        f_new[j_bot][i][7] = f8 + (f1 - f3)/2.0 - 0.5*rho_bot*U_bot + 1.0/6.0*rho_bot*V_bot;
    }
}