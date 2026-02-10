void bounce_back(
    int Nx,
    int Ny,
    int Q,
    double f[Ny][Nx][Q],
    double f_new[Ny][Nx][Q],
    int solid_mask[Ny][Nx],
    int cx[],
    int cy[],
    int opp[]
)
{
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {

            if (solid_mask[j][i]) continue;

            for (int k = 0; k < Q; k++){
                int in = i + (int)cx[k];
                int jn = j + (int)cy[k];

                if (in < 0 || in >= Nx || jn < 0 || jn >= Ny) continue;
                if (!solid_mask[jn][in]) continue;

                f_new[j][i][opp[k]] = f[j][i][k];
            }
        }
    }
}
