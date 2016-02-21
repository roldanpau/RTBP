#include <rtbp.h> 	// DIM

double re_integrand_B(double s, void *params);
int re_B(double mu, double p[DIM], double zu[DIM], double omega, double *res,
      int M, int N);
int im_B(double mu, double p[DIM], double zu[DIM], double omega, double *res,
      int M, int N);
int re_C(double mu, double p[DIM], double zs[DIM], double omega, double *res,
      int M, int N);
int im_C(double mu, double p[DIM], double zs[DIM], double omega, double *res,
      int M, int N);
