#include <rtbp.h>	// DIM
int orbitp_del(double mu, double p[DIM], int n, double *orbit);
int orbitp_del_bwd(double mu, double p[DIM], int n, double *orbit);
inline void print_orbit_del(int n, double *orbit);
