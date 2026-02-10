void interpolated_bounce_back(
    int Nx,
    int Ny,
    int Q,
    double f[Ny][Nx][Q],
    double f_new[Ny][Nx][Q],
    int solid_mask[Ny][Nx],
    int has_wall[Ny][Nx][Q],
    double q[Ny][Nx][Q],
    int opp[],
    int cx[],
    int cy[]
)
{
    for (int j = 0; j < Ny; j++) {
        for (int i = 0; i < Nx; i++) {

            for (int k = 0; k < Q; k++) {

                if (!has_wall[j][i][k]) continue;

                int is = i - cx[k];
                int js = j - cy[k];

                if (is < 0 || is >= Nx || js < 0 || js >= Ny) continue;

                double qq = q[j][i][k];

                if (qq <= 0.5) {
                    f_new[j][i][opp[k]] =
                        (1.0 - 2.0*qq) * f[j][i][k]
                      + 2.0*qq * f[js][is][k];
                } else {
                    if (solid_mask[js][is]) {
                        f_new[j][i][opp[k]] = f[j][i][k];
                    } else {
                    f_new[j][i][opp[k]] =
                        (1.0 / (2.0*qq)) * f[j][i][k]
                      + (2.0*qq - 1.0) / (2.0*qq) * f[js][is][opp[k]];
                      }
                }
            }
        }
    }
}
