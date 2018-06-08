#include <rtbp.h> 	// DIM

double re_integrand_B_stoch(double s, void *params);
int re_B_stoch(double mu, double p[DIM], double zu[DIM], double omega, 
        double *res, int M, int N);
int im_B_stoch(double mu, double p[DIM], double zu[DIM], double omega, 
        double *res, int M, int N);
int re_C_stoch(double mu, double p[DIM], double zs[DIM], double omega, 
        double *res, int M, int N);
int im_C_stoch(double mu, double p[DIM], double zs[DIM], double omega, 
        double *res, int M, int N);
