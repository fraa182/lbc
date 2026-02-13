void bounce_back(
    int Nx,
    int Ny,
    int Q,
    double f[Ny][Nx][Q],
    double f_new[Ny][Nx][Q],
    int solid_mask[Ny][Nx],
    int cx[],
    int cy[],
    int opp[],
    int isperiodic_x,
    int isperiodic_y
)
{
    #pragma omp parallel for
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {

            // Avoid solid mask
            if (solid_mask[j][i]) continue;

            for (int k = 0; k < Q; k++){
                int in = i + cx[k];
                int jn = j + cy[k];

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

                // Avoid solid mask
                if (!solid_mask[jn][in]) continue;

                // Apply bounce-back (BB)
                f_new[j][i][opp[k]] = f[j][i][k];
            }
        }
    }
}
