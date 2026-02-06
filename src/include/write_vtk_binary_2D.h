#ifndef WRITE_VTK_BINARY_2D_H
#define WRITE_VTK_BINARY_2D_H

void write_vtk_binary_2D(
    const char *filename,
    int Nx,
    int Ny,
    double dx,
    double u[Ny][Nx],
    double v[Ny][Nx],
    double rho[Ny][Nx]
);

#endif
