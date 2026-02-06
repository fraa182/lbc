#ifndef BOUNCE_BACK_H
#define BOUNCE_BACK_H

void bounce_back(
    int Nx,
    int Ny,
    int Q,
    double f_new[Ny][Nx][Q],
    int solid_mask[Ny][Nx],
    int opp[]
);

#endif
