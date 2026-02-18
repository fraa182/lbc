void streaming(
    int Nx,
    int Ny,
    int Q,
    double f[Ny][Nx][Q],
    double f_new[Ny][Nx][Q],
    int solid_mask[Ny][Nx],
    int cx[],
    int cy[],
    int isperiodic_x,
    int isperiodic_y
)
{
    #pragma omp parallel for
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {

            if (solid_mask[j][i]) continue;

            for (int k = 0; k < Q; k++) {
                int i_src = i - cx[k];
                int j_src = j - cy[k];

                // If we're using the periodic BC along x, wrap x
                if (isperiodic_x){
                    if (i_src < 0) i_src += Nx;
                    else if (i_src >= Nx) i_src -= Nx;
                }

                // If we're using the periodic BC along y, wrap y
                if (isperiodic_y){
                    if (j_src < 0) j_src += Ny;
                    else if (j_src >= Ny) j_src -= Ny;
                }

                // If non-periodic and out of bounds: leave it for BCs to set later
                if (i_src < 0 || i_src >= Nx || j_src < 0 || j_src >= Ny) continue;

                // Source is solid: this population must be provided by bounce-back/BC
                if (solid_mask[j_src][i_src]) continue;

                f_new[j][i][k] = f[j_src][i_src][k];
            }
        }
    }
}