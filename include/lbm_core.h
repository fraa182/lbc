#ifndef LBM_CORE_H
#define LBM_CORE_H

#include "lbm_types.h"

void main_lbm(
    int Nx,
    int Ny,
    int Q,
    double f[Ny][Nx][Q], 
    double f_new[Ny][Nx][Q], 
    double rho[Ny][Nx], 
    double u[Ny][Nx], 
    double v[Ny][Nx], 
    int cx[], 
    int cy[], 
    double w[], 
    int opp[], 
    double omega_eff[Ny][Nx], 
    int solid_mask[Ny][Nx], 
    Boundary boundaries[],
    int num_boundaries,
    int isperiodic_x,
    int isperiodic_y,
    double Fx,
    double Fy
);

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
    double w[],
    double Fx,
    double Fy
);

void streaming(
    int Nx,
    int Ny,
    int Q,
    double f[Ny][Nx][Q],
    double f_new[Ny][Nx][Q],
    int solid_mask[Ny][Nx],
    int cx[],
    int cy[]
);

void compute_macroscopic_fields(
    int Nx,
    int Ny,
    int Q,
    double f_new[Ny][Nx][Q],
    int solid_mask[Ny][Nx],
    int cx[],
    int cy[],
    double rho[Ny][Nx],
    double u[Ny][Nx],
    double v[Ny][Nx],
    double Fx,
    double Fy
);

#endif