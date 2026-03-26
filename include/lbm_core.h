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
    double Fy,
    double *L,
    double *D,
    double phi[Ny][Nx],
    int use_IBB
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
    int cy[],
    int isperiodic_x,
    int isperiodic_y
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

void compute_force_mem(
    int Nx, 
    int Ny,
    int Q,
    double f[Ny][Nx][Q],
    double f_new[Ny][Nx][Q],
    int solid_mask[Ny][Nx],
    int isperiodic_x,
    int isperiodic_y,
    const int cx[Q], 
    const int cy[Q],
    const int opp[Q],
    double *Fx, 
    double *Fy
);

void mei_initialization(
    int Nx,
    int Ny, 
    int Q,
    double (**pf)[Nx][Q], 
    double (**pf_new)[Nx][Q], 
    double rho_in[Ny][Nx],
    double U_in[Ny][Nx], 
    double V_in[Ny][Nx],
    double rho[Ny][Nx], 
    double u[Ny][Nx], 
    double v[Ny][Nx], 
    int cx[], 
    int cy[], 
    double w[], 
    double omega_eff[Ny][Nx], 
    int isperiodic_x,
    int isperiodic_y,
    int init_iter_max,
    double tol_rho
);

void sponge_layer(
    int Nx,
    int Ny,
    double dx,
    double xc,
    double yc,
    double r_inner,
    double r_outer,
    double tau_base,
    double incr_tau,
    double exp_coeff,
    double omega_eff[Ny][Nx]
);

#endif