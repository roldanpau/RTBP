#include <prtbp.h>	// section_t

extern const int ERR_DPRTBP_2D;
extern const int ERR_NOTREAL;

int hyper(double mu, section_t sec, double H, int n, double p[2], double eval[2], 
      double evec[4]);
