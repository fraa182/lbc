void periodic_y(
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
    
    for (int k = 0; k < Q; k++) {

        if (cy[k] > 0) {
            int j = 0;
            int j_src = Ny - cy[k];

            for (int i = 0; i < Nx; i++) {
                int i_src = i - cx[k];

                if (i_src < 0 || i_src >= Nx) continue;
                if (solid_mask[j][i]) continue;
                if (solid_mask[j_src][i_src]) continue;

                f_new[j][i][k] = f[j_src][i_src][k];
            }
        }

        if (cy[k] < 0) {
            int j = Ny - 1;
            int j_src = -cy[k] - 1;

            for (int i = 0; i < Nx; i++) {
                int i_src = i - cx[k];

                if (i_src < 0 || i_src >= Nx) continue;
                if (solid_mask[j][i]) continue;
                if (solid_mask[j_src][i_src]) continue;

                f_new[j][i][k] = f[j_src][i_src][k];
            }
        }
    }
}