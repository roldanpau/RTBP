#define DIM 4	// dimension of the (planar) RTBP
#define ERR_COLLISION 1
double Hamilt(double mu, const double *p);
int rtbp(double t, const double *x, double *y, void *params);
int rtbp_inv(double t, const double *x, double *y, void *params);
