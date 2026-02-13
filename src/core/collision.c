void collision(
    int Nx,
    int Ny,
    int Q,
    double f[Ny][Nx][Q],
    double rho[Ny][Nx],
    double u[Ny][Nx],
    double v[Ny][Nx],
    double omega_eff[Ny][Nx],
    int solid_mask[Ny][Nx],
    int cx[],
    int cy[],
    double w[],
    double Fx,
    double Fy
)
{
    #pragma omp parallel for
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {

            if (solid_mask[j][i]) continue;

            double ux = u[j][i];
            double uy = v[j][i];

            double v_sq = ux*ux + uy*uy;
            double uF = ux*Fx + uy*Fy;

            for (int k = 0; k < Q; k++) {
                double cu  = cx[k]*ux + cy[k]*uy;
                double feq = w[k]*rho[j][i]*(1.0 + 3.0*cu + 4.5*cu*cu - 1.5*v_sq);

                double cF = cx[k]*Fx + cy[k]*Fy;
                double Fk = w[k] * (3.0*cF + 9.0*cu*cF - 3.0*uF);

                f[j][i][k] += -omega_eff[j][i] * (f[j][i][k] - feq) + (1 - 0.5*omega_eff[j][i])*Fk;
            }
        }
    }
}
