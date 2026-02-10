#ifndef LBM_TYPES_H
#define LBM_TYPES_H

typedef void (*BCFunction)(
    int Nx, 
    int Ny, 
    int Q, 
    int index, 
    double val1, double val2,
    double f_new[Ny][Nx][Q],
    double rho[Ny][Nx], 
    double u[Ny][Nx], 
    double v[Ny][Nx], 
    int solid_mask[Ny][Nx], 
    int cx[],
    int cy[],
    double w[]
);

typedef struct {
    int index;
    double val1;
    double val2;
    BCFunction apply; 
} Boundary;

#endif