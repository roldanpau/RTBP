#include <rtbp.h>	// DIM
#include <frtbp.h>	// DIMV
#include <prtbp.h>	// section_t

extern const int ERR_VECTORFIELD;
extern const int ERR_DPRTBP;
extern const int ERR_POINCARE;

int set_dprtbp_2d(double dp[DIMV], double f[DIM], double g[DIM], double
      dp2d[4]);

int dprtbp_2d(double mu, section_t sec, double H, int cuts, double p[2], double dp2d[4]);
int dprtbp_2d_inv(double mu, section_t sec, double H, int cuts, double p[2], 
      double dp2d[4]);
