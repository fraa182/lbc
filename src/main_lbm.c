#include <stdio.h>
#include <math.h>
#include "include/lbm.h"
#define Q 9

void main_lbm(
    int Nx,
    int Ny, 
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
    int i_in, 
    double U_in, 
    int i_out, 
    double rho_out)
{

    // Collision step
    collision(Nx,Ny,Q,f,rho,u,v,omega_eff,solid_mask,cx,cy,w);

    // Bounce-back boundary condition at solid walls
    bounce_back(Nx,Ny,Q,f_new,solid_mask,opp);

    // Streaming step
    streaming(Nx,Ny,Q,f,f_new,cx,cy);

    // Velocity inlet
    velocity_inlet(Nx,Ny,Q,i_in,U_in,rho_out,f_new,u,v,solid_mask,cx,cy,w);

    // Pressure outlet
    pressure_outlet(Nx,Ny,Q,i_out,rho_out,f_new,rho,u,v,solid_mask,cx,cy,w);
    
    // Compute macroscopic quantities
    compute_macroscopic_fields(Nx,Ny,Q,f_new,solid_mask,cx,cy,rho,u,v);

}
