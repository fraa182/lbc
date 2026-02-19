#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "lbm_utils.h"

static inline uint32_t bswap32_u32(uint32_t x) {
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_bswap32(x);
#else
    return (x>>24) | ((x>>8)&0x0000FF00u) | ((x<<8)&0x00FF0000u) | (x<<24);
#endif
}

static inline float swap_float(float x) {
#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
    uint32_t u;
    memcpy(&u, &x, sizeof u);
    u = bswap32_u32(u);
    memcpy(&x, &u, sizeof x);
#endif
    return x;
}

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
)    
{
    FILE *fd = fopen(filename, "wb");
    if (!fd) {
        perror("Failed to open VTK file");
        return;
    }

    static char iobuf[1<<20];
    setvbuf(fd, iobuf, _IOFBF, sizeof(iobuf));

    /* ---- Header (ASCII) ---- */
    fprintf(fd, "# vtk DataFile Version 3.0\n");
    fprintf(fd, "2D LBM output\n");
    fprintf(fd, "BINARY\n");
    fprintf(fd, "DATASET STRUCTURED_POINTS\n");
    fprintf(fd, "DIMENSIONS %d %d 1\n", Nx, Ny);
    fprintf(fd, "ORIGIN 0 0 0\n");
    fprintf(fd, "SPACING %f %f 1\n", dx, dx);
    fprintf(fd, "\nPOINT_DATA %d\n", Nx * Ny);

    /* ---- Allocate velocity and density ---- */
    float (*velocity_write)[Nx][3] = malloc(Ny * sizeof *velocity_write);
    float (*rho_write)[Nx] = malloc(Ny * sizeof *rho_write);

    /* ---- Fill velocity and density ---- */
    #pragma omp parallel for
    for (int j = 0; j < Ny; j++) {
        for (int i = 0; i < Nx; i++) {
            velocity_write[j][i][0] = swap_float(u[j][i]*cu);
            velocity_write[j][i][1] = swap_float(v[j][i]*cu);
            velocity_write[j][i][2] = swap_float(0.0f);
            rho_write[j][i] = swap_float(rho[j][i]*crho);
        }
    }

    /* ---- Write Velocity field ---- */
    fprintf(fd, "VECTORS velocity float\n");
    fwrite(velocity_write, sizeof(float), Nx*Ny*3, fd);

    /* ---- Write Density field ---- */
    fprintf(fd, "\nSCALARS density float 1\n");
    fprintf(fd, "LOOKUP_TABLE default\n");
    fwrite(rho_write, sizeof(float), Nx*Ny, fd);

    free(velocity_write);
    free(rho_write);

    fclose(fd);
}
