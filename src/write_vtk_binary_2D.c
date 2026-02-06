#include <stdio.h>
#include "include/swap_float.h"

void write_vtk_binary_2D(
    const char *filename,
    int Nx,
    int Ny,
    double dx,
    double u[Ny][Nx],
    double v[Ny][Nx],
    double rho[Ny][Nx]
)    
{
    FILE *fd = fopen(filename, "wb");
    if (!fd) {
        perror("Failed to open VTK file");
        return;
    }

    /* ---- Header (ASCII) ---- */
    fprintf(fd, "# vtk DataFile Version 3.0\n");
    fprintf(fd, "2D LBM output\n");
    fprintf(fd, "BINARY\n");
    fprintf(fd, "DATASET STRUCTURED_POINTS\n");
    fprintf(fd, "DIMENSIONS %d %d 1\n", Nx, Ny);
    fprintf(fd, "ORIGIN 0 0 0\n");
    fprintf(fd, "SPACING %f %f 1\n", dx, dx);
    fprintf(fd, "\nPOINT_DATA %d\n", Nx * Ny);

    /* ---- Velocity field ---- */
    fprintf(fd, "VECTORS velocity float\n");

    for (int j = 0; j < Ny; j++) {
        for (int i = 0; i < Nx; i++) {
            float vec[3];
            vec[0] = swap_float(u[j][i]);
            vec[1] = swap_float(v[j][i]);
            vec[2] = swap_float(0.0f);
            fwrite(vec, sizeof(float), 3, fd);
        }
    }

    /* ---- Density field ---- */
    fprintf(fd, "\nSCALARS density float 1\n");
    fprintf(fd, "LOOKUP_TABLE default\n");

    for (int j = 0; j < Ny; j++) {
        for (int i = 0; i < Nx; i++) {
            float r = swap_float(rho[j][i]);
            fwrite(&r, sizeof(float), 1, fd);
        }
    }

    fclose(fd);
}
