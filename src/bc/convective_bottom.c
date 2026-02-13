void convective_bottom(
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
    int j_bottom = index;

    for (int i = 0; i < Nx; i++) {

        if (solid_mask[j_bottom][i]) continue;

        for (int k = 0; k < Q; k++){
            f_new[j_bottom][i][k] = f_new[j_bottom+1][i][k];
        }
    }
}