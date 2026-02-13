void velocity_top(
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
    double U_top = val1;
    double V_top = val2;
    int j_top = index;
    
    #pragma omp parallel for
    for (int i = 0; i < Nx; i++) {

        if (solid_mask[j_top][i]) continue;

        double f0 = f_new[j_top][i][0];
        double f1 = f_new[j_top][i][1];
        double f2 = f_new[j_top][i][2];
        double f3 = f_new[j_top][i][3]; 
        double f5 = f_new[j_top][i][5];
        double f6 = f_new[j_top][i][6];

        double rho_top = (f0 + f1 + f3 + 2*(f2 + f5 + f6))/(1 + V_top);

        f_new[j_top][i][4] = f2 - 2.0/3.0*rho_top*V_top;
        f_new[j_top][i][7] = f5 + (f1 - f3)/2.0 - 0.5*rho_top*U_top - 1/6.0*rho_top*V_top;
        f_new[j_top][i][8] = f6 - (f1 - f3)/2.0 + 0.5*rho_top*U_top - 1/6.0*rho_top*V_top;

        rho[j_top][i] = rho_top;
        u[j_top][i]   = U_top;
        v[j_top][i]   = V_top;
    }
}