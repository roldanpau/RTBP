#include <rtbp.h>	// DIM

int intersec_unst(double mu, double H, double p[2], double v[2], 
      double lambda, int n, double h1, double h2, double l, 
      double *h, double p_u[DIM], double *t, double z[DIM]);
int intersec_st(double mu, double H, double p[2], double v[2], 
      double lambda, int n, double h1, double h2, double l,
      double *h, double p_s[DIM], double *t, double z[DIM]);
