#ifndef VELOCITY_INLET_H
#define VELOCITY_INLET_H

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
);

#endif
