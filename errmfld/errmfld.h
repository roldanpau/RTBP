#include <prtbp.h>	// section_t

double err_mfld(double mu, section_t sec, double H, int k, double p[2],
      double v[2], double lambda, int stable, double h);
double h_opt(double mu, section_t sec, double H, int k, double p[2], 
      double v[2], double lambda, int stable);
