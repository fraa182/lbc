#ifndef COLLISION_H
#define COLLISION_H

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
    double w[]
);

#endif
