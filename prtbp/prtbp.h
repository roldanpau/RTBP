#ifndef PRTBP_H_INCLUDED
#define PRTBP_H_INCLUDED

#include <section.h>	// section_t
#include "rtbp.h"	// DIM

extern const double POINCARE_TOL;	///< error bound (tolerance) for Poincare map

int prtbp(double mu, section_t sec, int cuts, double x[DIM], double *ti);
int prtbp_inv(double mu, section_t sec, int cuts, double x[DIM], double *ti);

#endif // PRTBP_H_INCLUDED
