#ifndef STREAMING_H
#define STREAMING_H

void streaming(
    int Nx,
    int Ny,
    int Q,
    double f[Ny][Nx][Q],
    double f_new[Ny][Nx][Q],
    int cx[],
    int cy[]
);

#endif
