void bounce_back(
    int Nx,
    int Ny,
    int Q,
    double f_new[Ny][Nx][Q],
    int solid_mask[Ny][Nx],
    int opp[]
)
{
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {

            if (!solid_mask[j][i]) continue;

            double tmp[Q];
            for (int k = 0; k < Q; k++)
                tmp[k] = f_new[j][i][k];

            for (int k = 0; k < Q; k++)
                f_new[j][i][k] = tmp[opp[k]];
        }
    }
}
