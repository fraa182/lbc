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
    for (int k = 0; k < Q; k++) {
        for (int i = 0; i < Nx; i++) {

            int i_src = i - (int)cx[k];
            if (i_src < 0)      i_src += Nx;
            if (i_src >= Nx)    i_src -= Nx;

            for (int j = 0; j < Ny; j++) {

                int j_src = j - (int)cy[k];
                if (j_src < 0)      j_src += Ny;
                if (j_src >= Ny)    j_src -= Ny;

                if (solid_mask[j_src][i_src]) continue;
                if (solid_mask[j][i]) continue;

                f_new[j][i][k] = f[j_src][i_src][k];
            }
        }
    }
}