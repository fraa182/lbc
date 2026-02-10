#ifndef LBM_UTILS_H
#define LBM_UTILS_H

int ensure_directory_exists(const char *path);

float swap_float(float val);

void write_vtk_binary_2D(
    const char *filename,
    int Nx,
    int Ny,
    double dx,
    double u[Ny][Nx],
    double v[Ny][Nx],
    double rho[Ny][Nx],
    double cu,
    double crho
);

#endif