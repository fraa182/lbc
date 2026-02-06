#ifndef PRESSURE_OUTLET_H
#define PRESSURE_OUTLET_H

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
);

#endif
