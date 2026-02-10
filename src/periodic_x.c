void periodic_x(
    int Nx, 
    int Ny, 
    int Q, 
    int index, 
    double v1, 
    double v2,
    double f_new[Ny][Nx][Q], 
    double rho[Ny][Nx], 
    double u[Ny][Nx], 
    double v[Ny][Nx], 
    int solid_mask[Ny][Nx], 
    int cx[], 
    int cy[], 
    double w[]
) {
    
    for (int j = 0; j < Ny; j++) {

        for (int k = 0; k < Q; k++) {
            if (cx[k] < 0) {
                f_new[j][0][k] = f_new[j][Nx-1][k];
            }
            if (cx[k] > 0) {
                f_new[j][Nx-1][k] = f_new[j][0][k];
            }
        }
    }
}