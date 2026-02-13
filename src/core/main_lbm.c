#include <stdio.h>
#include <math.h>
#include "lbm.h"

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
    double Fy)
{

    // Collision step
    collision(Nx,Ny,Q,f,rho,u,v,omega_eff,solid_mask,cx,cy,w,Fx,Fy);

    // Bounce-back boundary condition at solid walls
    bounce_back(Nx,Ny,Q,f,f_new,solid_mask,cx,cy,opp,isperiodic_x,isperiodic_y);

    // Streaming step
    streaming(Nx,Ny,Q,f,f_new,solid_mask,cx,cy);

    // Boundary conditions
    for (int i = 0; i < num_boundaries; i++) {
        boundaries[i].apply(
            Nx, Ny, Q, 
            boundaries[i].index, 
            boundaries[i].val1, boundaries[i].val2,
            f, f_new, rho, u, v, solid_mask, cx, cy, w
        );
    }
    
    // Compute macroscopic quantities
    compute_macroscopic_fields(Nx,Ny,Q,f_new,solid_mask,cx,cy,rho,u,v,Fx,Fy);

}
