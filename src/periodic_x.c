void periodic_x(
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

        if (cx[k] > 0) {
            int i = 0;
            int i_src = Nx - cx[k];

            for (int j = 0; j < Ny; j++) {
                int j_src = j - cy[k];

                if (j_src < 0 || j_src >= Ny) continue;
                if (solid_mask[j][i]) continue;
                if (solid_mask[j_src][i_src]) continue;

                f_new[j][i][k] = f[j_src][i_src][k];
            }
        }

        if (cx[k] < 0) {
            int i = Nx - 1;
            int i_src = -cx[k] - 1;

            for (int j = 0; j < Ny; j++) {
                int j_src = j - cy[k];

                if (j_src < 0 || j_src >= Ny) continue;
                if (solid_mask[j][i]) continue;
                if (solid_mask[j_src][i_src]) continue;

                f_new[j][i][k] = f[j_src][i_src][k];
            }
        }
    }
}