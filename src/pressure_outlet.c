#define myzero 1e-16

void pressure_outlet(
    int Nx,
    int Ny,
    int Q,
    int i_out,
    double rho_out,
    double f_new[Ny][Nx][Q],
    double rho[Ny][Nx],
    double u[Ny][Nx],
    double v[Ny][Nx],
    int solid_mask[Ny][Nx],
    int cx[],
    int cy[],
    double w[]
)
{
    int i_int = (i_out > 0) ? i_out - 1 : i_out + 1;

    for (int j = 0; j < Ny; j++) {

        if (solid_mask[j][i_out]) continue;

        double s  = 0.0;
        double sx = 0.0;
        double sy = 0.0;

        for (int k = 0; k < Q; k++) {
            double fk = f_new[j][i_int][k];
            s  += fk;
            sx += fk * cx[k];
            sy += fk * cy[k];
        }

        if (s < myzero) s = myzero;

        double ux = sx / s;
        double uy = sy / s;
        double v_sq = ux*ux + uy*uy;

        rho[j][i_out] = rho_out;
        u[j][i_out]   = ux;
        v[j][i_out]   = uy;

        for (int k = 0; k < Q; k++) {
            double cu = cx[k]*ux + cy[k]*uy;
            f_new[j][i_out][k] =
                w[k]*rho_out *
                (1.0 + 3.0*cu + 4.5*cu*cu - 1.5*v_sq);
        }
    }
}
