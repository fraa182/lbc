void convective_outlet(
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
    int i_out = index;

    #pragma omp parallel for
    for (int j = 0; j < Ny; j++) {

        if (solid_mask[j][i_out]) continue;

        for (int k = 0; k < Q; k++){
            f_new[j][i_out][k] = f_new[j][i_out-1][k];
        }
    }
}
