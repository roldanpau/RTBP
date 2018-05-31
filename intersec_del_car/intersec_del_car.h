int intersec_del_car_unst(double mu, section_t sec, branch_t br, double H,
        double p[2], double v[2], double lambda, int n, double h1, double h2,
        double l, double *h, double p_u[2], double *t, double z_del[DIM], 
        double z_car[DIM], double z_u[DIM], double z_u_car[DIM]);
int intersec_del_car_st(double mu, section_t sec, branch_t br, double H,
        double p[2], double v[2], double lambda, int n, double h1, double h2,
        double l, double *h, double p_u[2], double *t, double z_del[DIM], 
        double z_car[DIM], double z_u[DIM], double z_u_car[DIM]);
