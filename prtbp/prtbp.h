#ifndef PRTBP_H_INCLUDED
#define PRTBP_H_INCLUDED

#include "rtbp.h"	// DIM

// Different Poincare sections
// SEC1 is {y=0, p_y>0}, SEC2 is {y=0, p_y<0}
typedef enum {SEC1, SEC2} section_t;    

int prtbp(double mu, section_t sec, int cuts, double x[DIM], double *ti);
int prtbp_inv(double mu, section_t sec, int cuts, double x[DIM], double *ti);

#endif // PRTBP_H_INCLUDED
