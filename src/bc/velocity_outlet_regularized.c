void velocity_outlet_regularized(
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
    int i_out = index;
    int i_inter = index - 1;

    #pragma omp parallel for
    for (int j = 0; j < Ny; j++) {

        if (solid_mask[j][i_out]) continue;

        // Retrieve interior quantities
        double rho_inter = rho[j][i_inter];
        double u_inter = u[j][i_inter];
        double v_inter = v[j][i_inter];

        double v_sq_inter = u_inter*u_inter + v_inter*v_inter;

        // Retrieve outlet quantities
        double rho_out = rho_inter;
        double U_out = val1[j]; 
        double V_out = val2[j];

        double v_sq_out = U_out*U_out + V_out*V_out;

        // Extrapolate stress tensor from interior nodes
        double Pi_xx = 0.0;
        double Pi_xy = 0.0;
        double Pi_yy = 0.0;
        
        for (int k = 0; k < Q; k++) {
            double cu = cx[k]*u_inter + cy[k]*v_inter;
            double feq = w[k] * rho_inter * (1.0 + 3.0*cu + 4.5*cu*cu - 1.5*v_sq_inter);
            double fneq = f_new[j][i_inter][k] - feq;

            Pi_xx += cx[k]*cx[k]*fneq;
            Pi_xy += cx[k]*cy[k]*fneq;
            Pi_yy += cy[k]*cy[k]*fneq;
        }

        // Update populations
        for (int k = 0; k < Q; k++) {
            // Compute equilbirium
            double cu_out = cx[k]*U_out + cy[k]*V_out;
            double feq_out = w[k] * rho_out * (1.0 + 3.0*cu_out + 4.5*cu_out*cu_out - 1.5*v_sq_out);

            // Compute regularized non-equilibrium with 2nd order Hermite polynomials
            double Qxx = cx[k]*cx[k] - 1.0/3.0;
            double Qxy = cx[k]*cy[k];
            double Qyy = cy[k]*cy[k] - 1.0/3.0;

            double fneq_reg = (3.0/2.0) * w[k] * (Qxx*Pi_xx + 2.0*Qxy*Pi_xy + Qyy*Pi_yy);

            // Sum equilibrium and non-equilibrium populations
            f_new[j][i_out][k] = feq_out + fneq_reg;
        }
    }
}