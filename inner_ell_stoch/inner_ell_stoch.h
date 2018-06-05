double re_integrand_inner_ell_stoch(double s, void *params);
double re_f_integrand_stoch(double mu, double x2[DIM]);
double im_f_integrand_stoch(double mu, double x2[DIM]);
int re_inner_ell_stoch(double mu, double T, double x[DIM], double *re_A);
int im_inner_ell_stoch(double mu, double T, double x[DIM], double *im_A);
