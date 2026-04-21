void pressure_top_regularized(
    int Nx, 
    int Ny,
    int Q, 
    int index,
    const double *val1, 
    const double *val2,
    double f[Ny][Nx][Q],
    double f_new[Ny][Nx][Q], 
    double rho[Ny][Nx], 
    double u[Ny][Nx], 
    double v[Ny][Nx], 
    int solid_mask[Ny][Nx], 
    int cx[],
    int cy[],
    double w[]
){
    int j_top = index;
    int j_inter = index - 1;

    #pragma omp parallel for
    for (int i = 0; i < Nx; i++) {
        if (solid_mask[j_top][i]) continue;

        // Retrieve interior quantities
        double rho_inter = rho[j_inter][i];
        double u_inter = u[j_inter][i];
        double v_inter = v[j_inter][i];

        double v_sq_inter = u_inter*u_inter + v_inter*v_inter;

        // Retrieve top quantities
        double rho_top = val1[i];
        double U_top = u_inter; 
        double V_top = v_inter;

        double v_sq_top = U_top*U_top + V_top*V_top;

        // Extrapolate stress tensor from interior nodes
        double Pi_xx = 0.0;
        double Pi_xy = 0.0;
        double Pi_yy = 0.0;
        
        for (int k = 0; k < Q; k++) {
            double cu = cx[k]*u_inter + cy[k]*v_inter;
            double feq = w[k] * rho_inter * (1.0 + 3.0*cu + 4.5*cu*cu - 1.5*v_sq_inter);
            double fneq = f_new[j_inter][i][k] - feq;

            Pi_xx += cx[k]*cx[k]*fneq;
            Pi_xy += cx[k]*cy[k]*fneq;
            Pi_yy += cy[k]*cy[k]*fneq;
        }

        // Update populations
        for (int k = 0; k < Q; k++) {
            // Compute equilbirium
            double cu_top = cx[k]*U_top + cy[k]*V_top;
            double feq_top = w[k] * rho_top * (1.0 + 3.0*cu_top + 4.5*cu_top*cu_top - 1.5*v_sq_top);

            // Compute regularized non-equilibrium with 2nd order Hermite polynomials
            double Qxx = cx[k]*cx[k] - 1.0/3.0;
            double Qxy = cx[k]*cy[k];
            double Qyy = cy[k]*cy[k] - 1.0/3.0;

            double fneq_reg = (3.0/2.0) * w[k] * (Qxx*Pi_xx + 2.0*Qxy*Pi_xy + Qyy*Pi_yy);

            // Sum equilibrium and non-equilibrium populations
            f_new[j_top][i][k] = feq_top + fneq_reg;
        }
    }
}