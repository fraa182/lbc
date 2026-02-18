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
){
    double fx = 0.0;
    double fy = 0.0;

    #pragma omp parallel for reduction(+:fx,fy)
    for (int j = 0; j < Ny; ++j){
        for (int i = 0; i < Nx; ++i){

            if (solid_mask[j][i]) continue; // fluid nodes only

            for (int k = 1; k < Q; ++k){ // skip k=0 (rest)
                int in = i + cx[k];
                int jn = j + cy[k];

                // If we're using the periodic BC along x, wrap x
                if (isperiodic_x) {
                    if (in < 0) in += Nx;
                    else if (in >= Nx) in -= Nx;
                } else {
                    if (in < 0 || in >= Nx) continue;
                }

                // If we're using the periodic BC along y, wrap y
                if (isperiodic_y){
                    if (jn < 0) jn += Ny;
                    else if (jn >= Ny) jn -= Ny;
                } else {
                    if (jn < 0 || jn >= Ny) continue;
                }

                if (solid_mask[jn][in]) {
                    int ko = opp[k];

                    // Momentum exchange across the fluid-solid link
                    double df = f[j][i][k] + f_new[j][i][ko];
                    fx += df * cx[k];
                    fy += df * cy[k];
                }
            }
        }
    }

    *Fx = fx;
    *Fy = fy;
}
