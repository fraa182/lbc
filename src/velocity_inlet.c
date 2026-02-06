#define myzero 1e-16

void velocity_inlet(
    int Nx,
    int Ny,
    int Q,
    int i_in,
    double U_in,
    double rho_out,
    double f_new[Ny][Nx][Q],
    double u[Ny][Nx],
    double v[Ny][Nx],
    int solid_mask[Ny][Nx],
    int cx[],
    int cy[],
    double w[]
)
{
    for (int j = 0; j < Ny; j++) {

        if (solid_mask[j][i_in]) continue;

        u[j][i_in] = U_in;
        v[j][i_in] = 0.0;

        double v_sq = U_in * U_in;

        for (int k = 0; k < Q; k++) {
            double cu = cx[k]*U_in;
            f_new[j][i_in][k] =
                w[k]*rho_out *
                (1.0 + 3.0*cu + 4.5*cu*cu - 1.5*v_sq);
        }
    }
}
