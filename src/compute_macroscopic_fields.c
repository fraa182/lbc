#include <math.h>
#define myzero 1e-16

void compute_macroscopic_fields(
    int Nx,
    int Ny,
    int Q,
    double f_new[Ny][Nx][Q],
    int solid_mask[Ny][Nx],
    int cx[],
    int cy[],
    double rho[Ny][Nx],
    double u[Ny][Nx],
    double v[Ny][Nx],
    double Fx,
    double Fy
)
{
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {

            double s  = 0.0;
            double sx = 0.0;
            double sy = 0.0;

            for (int k = 0; k < Q; k++) {
                double fk = f_new[j][i][k];
                s  += fk;
                sx += fk * cx[k];
                sy += fk * cy[k];
            }

            if (s < myzero) {
                s = myzero;
            }

            if (!solid_mask[j][i]) {
                rho[j][i] = s;
                u[j][i]   = (sx + 0.5*Fx) / s;
                v[j][i]   = (sy + 0.5*Fy) / s;
            } else {
                rho[j][i] = NAN;
                u[j][i]   = NAN;
                v[j][i]   = NAN;
            }
        }
    }
}
