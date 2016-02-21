#include <rtbp.h>	// DIM
#include <frtbp.h>	// DIMV
#include <prtbp.h>	// section_t
int dprtbp(double mu, section_t sec, int cuts, double x[DIM], double dp[DIMV]);
int dprtbp_inv(double mu, section_t sec, int cuts, double x[DIM], double dp[DIMV]);
