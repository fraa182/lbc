#include <stdio.h>
#include <math.h>
#include "../src/lbm_step.h"
#define Q 9

int main(){

    // Physical parameters
    int res = 10;
    double R = 0.02;
    double vx_size = R/res;

    double U_phys = 10.0;
    double Lx = 0.5;
    double Ly = 0.2;
    double t_final = 1;

    // Lattice domain
    int Nx = Lx/vx_size;
    int Ny = Ly/vx_size;

    double U_in = 0.1;
    double rho0 = 1.0;
    int i_in = 0;
    int i_out = Nx - 1;

    double cu = U_phys/U_in;
    double ct = vx_size/cu;

    double nu_phys = 0.000015*20;
    double nu_lattice = nu_phys/(cu*vx_size);
    double tau_base = 0.5 + 3*nu_lattice;

    int Nt = t_final/ct;

    int cx[Q] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
    int cy[Q] = {0, 0, 1, 0, -1, 1, 1, -1, -1};
    double w[Q] = {4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};
    int opp[Q] = {0, 3, 4, 1, 2, 7, 8, 5, 6};

    // Solid mask
    int solid_mask[Ny][Nx];
    for (int i = 0; i < Nx; i++){
        double x = i*vx_size;
        for (int j = 0; j < Ny; j++){
            double y = j*vx_size;
            solid_mask[j][i] = ((x - Lx/2)*(x - Lx/2) + (y - Ly/2)*(y - Ly/2)) <= R*R ? 1 : 0;
        }
    }

    // Collision operator
    double omega_eff[Ny][Nx];
    for (int i = 0; i < Nx; i++){
        for (int j = 0; j < Ny; j++){
            omega_eff[j][i] = 1.0/tau_base;
        }
    }

    // Flow field
    double rho[Ny][Nx], u[Ny][Nx], v[Ny][Nx];
    for (int i = 0; i < Nx; i++){
        for (int j = 0; j < Ny; j++){
            u[j][i] = 0;
            v[j][i] = 0;
            rho[j][i] = 1.0;
        }
    }
    
    // Particle distribution function
    double f[Ny][Nx][Q], f_new[Ny][Nx][Q];
    for (int i = 0; i < Nx; i++){
        for (int j = 0; j < Ny; j++){
            for (int k = 0; k < Q; k++){
                f[j][i][k] = w[k]*rho0;
                f_new[j][i][k] = 0;
            }
        }
    }

    // Temporal loop
    for (int it = 0; it < Nt; it++){
        // Print the timestep
        printf("Step %d of %d\n",it+1,Nt);

        // Execute the stream and collide LBM phases
        lbm_step(Nx,Ny,f,f_new,rho,u,v,cx,cy,w,opp,omega_eff,solid_mask,i_in,U_in,i_out,rho0);

        // Swap f with f_new
        for (int i = 0; i < Nx; i++){
            for (int j = 0; j < Ny; j++){
                for (int k = 0; k < Q; k++) {
                    double tmp = f[j][i][k];
                    f[j][i][k] = f_new[j][i][k];
                    f_new[j][i][k] = tmp;
                }
            }
        }

        // Write solution
        if (it % 100 == 0){
            char filename[256];
            snprintf(filename, sizeof(filename),"sol/fields_%05d.vtk", it);
            FILE *fd = fopen(filename, "w");

            fprintf(fd, "# vtk DataFile Version 3.0\n");
            fprintf(fd, "2D Grid Data\n");
            fprintf(fd, "ASCII\n");
            fprintf(fd, "DATASET STRUCTURED_GRID\n");
            fprintf(fd, "DIMENSIONS %d %d 1\n", Nx, Ny);
            fprintf(fd, "POINTS %d float\n", Nx*Ny);

            for (int j = 0; j < Ny; j++) {
                double y = j*vx_size;
                for (int i = 0; i < Nx; i++) {
                    double x = i*vx_size;
                    fprintf(fd, "%f %f 0.0\n", x, y);
                }
            }

            fprintf(fd, "\nPOINT_DATA %d\n", Nx*Ny);
            fprintf(fd, "VECTORS velocity float\n");
            for (int j=0; j < Ny; j++){
                for (int i=0; i < Nx; i++){
                    fprintf(fd, "%f %f 0.0\n", u[j][i], v[j][i]);
                }
            }

            fprintf(fd, "\nSCALARS density float 1\n");
            fprintf(fd, "LOOKUP_TABLE default\n");
            for (int j = 0; j < Ny; j++){
                for (int i = 0; i < Nx; i++){
                    fprintf(fd, "%f\n", rho[j][i]);
                }
            }

            fclose(fd);
        }

    }

    return 0;

}
