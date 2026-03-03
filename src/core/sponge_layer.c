#include <math.h>

void sponge_layer(
    int Nx,
    int Ny,
    double dx,
    double xc,
    double yc,
    double r_inner,
    double r_outer,
    double tau_base,
    double incr_tau,
    double exp_coeff,
    double omega_eff[Ny][Nx]
){

    // Initialize effective relaxation time
    double tau_eff = 0.0;

    // Compute omega, i.e., the inverse of the effective relaxation time, for each cell
    for (int j = 0; j < Ny; j++){
        double y = j*dx - yc;
        for (int i = 0; i < Nx; i++){
            double x = i*dx - xc;

            // Compute the position in the domain
            double r = sqrt(x*x + y*y);

            // Compute the effective relaxation time
            if (r <= r_inner){
                // Inside the inner (base) layer (no viscosity increment)
                tau_eff = tau_base;
            } else if ((r > r_inner) & (r <= r_outer)) {
                // Between inner and outer layer (exponential viscosity increment)
                tau_eff = tau_base + (incr_tau - 1)*tau_base * (exp(exp_coeff * (r - r_inner)) - 1)/(exp(exp_coeff * (r_outer - r_inner)) - 1);
            } else {
                // Outside the outer layer (incremented viscosity)
                tau_eff = tau_base*incr_tau;
            }

            // Update omega
            omega_eff[j][i] = 1.0/tau_eff;
        }
    }

}