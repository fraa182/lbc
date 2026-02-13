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
    int cy[],
    int isperiodic_x,
    int isperiodic_y
)
{
    #pragma omp parallel for
    for (int j = 0; j < Ny; j++) {
        for (int i = 0; i < Nx; i++) {

            // Avoid solid mask
            if (solid_mask[j][i]) continue;

            for (int k = 0; k < Q; k++) {

                if (!has_wall[j][i][k]) continue;

                int in = i - cx[k];
                int jn = j - cy[k];

                // If we're using the periodic BC along x, wrap x
                if (isperiodic_x) {
                    if (in < 0) in += Nx;
                    else if (in >= Nx) in -= Nx;
                } else {
                    if (in < 0 || in >= Nx) continue;
                }

                // If we're using the periodic BC along y, wrap y
                if (isperiodic_y){
                    if (jn < 0) jn += Ny;
                    else if (jn >= Ny) jn -= Ny;
                } else {
                    if (jn < 0 || jn >= Ny) continue;
                }

                // Apply IBB (Bouzidi)
                if (q[j][i][k] <= 0.5) {
                    f_new[j][i][opp[k]] =
                        (1.0 - 2.0*q[j][i][k]) * f[j][i][k]
                      + 2.0*q[j][i][k] * f[jn][in][k];
                } else {
                    if (solid_mask[jn][in]) {
                        f_new[j][i][opp[k]] = f[j][i][k];
                    } else {
                    f_new[j][i][opp[k]] =
                        (1.0 / (2.0*q[j][i][k])) * f[j][i][k]
                      + (2.0*q[j][i][k] - 1.0) / (2.0*q[j][i][k]) * f[jn][in][opp[k]];
                      }
                }
            }
        }
    }
}
