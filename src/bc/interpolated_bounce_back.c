#include <math.h>
#define eps 1e-12

void interpolated_bounce_back(
    int Nx,
    int Ny,
    int Q,
    double f[Ny][Nx][Q],
    double f_new[Ny][Nx][Q],
    int solid_mask[Ny][Nx],
    double phi[Ny][Nx],
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

                // Forward neighbor along +ck
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

                // Wall link?
                if (!solid_mask[jn][in]) continue;

                // Compute q from phi (signed distance)
                double denom = phi[j][i] - phi[jn][in];
                if (fabs(denom) < eps) {
                    // fallback to mid-link
                    denom = (denom >= 0.0) ? eps : -eps;
                }
                double q = phi[j][i]/denom;

                // Clamp q to a safe range
                if (q < eps) q = eps;
                if (q > 1.0 - eps) q = 1.0 - eps;

                // Behind node (fluid side) at x_f - c_k
                int ib = i - cx[k];
                int jb = j - cy[k];

                // If we're using the periodic BC along x, wrap x
                if (isperiodic_x) {
                    if (ib < 0) ib += Nx;
                    else if (ib >= Nx) ib -= Nx;
                } else {
                    if (ib < 0 || ib >= Nx){
                        // No valid behind node -> revert to simple BB
                        f_new[j][i][opp[k]] = f[j][i][k]; 
                        continue;
                    }
                }

                // If we're using the periodic BC along y, wrap y
                if (isperiodic_y){
                    if (jb < 0) jb += Ny;
                    else if (jb >= Ny) jb -= Ny;
                } else {
                    if (jb < 0 || jb >= Ny){
                        // No valid behind node -> revert to simple BB
                        f_new[j][i][opp[k]] = f[j][i][k]; 
                        continue;
                    }
                }

                // Apply IBB (Bouzidi)
                if (q <= 0.5) {
                    f_new[j][i][opp[k]] = (1.0 - 2.0*q) * f[j][i][k] + 2.0*q * f[jb][ib][k];
                } else {
                    if (solid_mask[jb][ib]) {
                        f_new[j][i][opp[k]] = f[j][i][k];
                    } else {
                        f_new[j][i][opp[k]] = (1.0 / (2.0*q)) * f[j][i][k] + (2.0*q - 1.0) / (2.0*q) * f[jb][ib][opp[k]];
                    }
                }
            }
        }
    }
}
