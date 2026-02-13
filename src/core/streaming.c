void streaming(
    int Nx,
    int Ny,
    int Q,
    double f[Ny][Nx][Q],
    double f_new[Ny][Nx][Q],
    int solid_mask[Ny][Nx],
    int cx[],
    int cy[]
)
{
    #pragma omp parallel for
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < Q; k++) {

                int i_src = i - cx[k];
                int j_src = j - cy[k];

                if (i_src < 0 || i_src >= Nx || j_src < 0 || j_src >= Ny) continue;
                if (solid_mask[j_src][i_src]) continue;
                if (solid_mask[j][i]) continue;

                f_new[j][i][k] = f[j_src][i_src][k];
            }
        }
    }
}