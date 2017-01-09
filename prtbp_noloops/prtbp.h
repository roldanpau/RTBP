/*! \file
    \brief Poincare map of the Restricted Three Body Problem

    $Author: roldan $
    $Date: 2013-03-26 22:22:14 $
*/

#ifndef PRTBP_H_INCLUDED
#define PRTBP_H_INCLUDED

#include <section.h>	// section_t
#include "rtbp.h"	// DIM

extern const double POINCARE_TOL;	///< error bound (tolerance) for Poincare map

/**
  Poincare map of the Restricted Three Body Problem.

  Consider the RTBP in rotating coordinates.
  Let S be the Poincare section SEC1 or SEC2, corresponding to {y=0, v_y>0}
  or {y=0, v_y<0}.
  Let $x$ be a point in the section, and suppose that flow at $x$ is
  transversal to the section.
  Compute the n-th iterate of the Poincare map $P^n(x)$ of the RTBP.
  This procedure also computes, as a side product, the integration time to
  intersect the Poincare section "n" times.

  \param[in] mu mass parameter for the RTBP
  \param[in] sec type of Poincare section (sec = SEC1 or SEC2).

  \param[in] cuts 
  number of iterates of the Poincare map (cuts=n: n cuts with the Poincare
  section).

  \param[in,out] x 
  Initial point, 4 coordinates: (X, Y, P_X, P_Y). 
  On return of the this function, it holds the image point $P^n(x)$.

  \param[out] ti
  On return, it holds the integration time to intersect the Poincare
  section "n" times.

  \return
  Returns a non-zero error code to indicate an error and 0 to indicate
  success.
  If an integration error is encountered, the function returns a non-zero
  value.

  \remark
  We are only interested in orbits that go around the origin. For some values
  of the energy, the orbit has "contractible" loops, which we want to avoid.
  For this, we impose that TWO CONSECUTIVE ITERATES do NOT lie both to the
  right or to the left of the origin.

  \remark
  On successful return of this function, the point $x$ is exactly on the
  section, i.e. we set coordinate $y$ exactly equal to zero.
*/    

int prtbp(double mu, section_t sec, int cuts, double x[DIM], double *ti);

/**
  Inverse Poincare map of the Restricted Three Body Problem.

  Consider the RTBP in rotating coordinates.
  Let S be the Poincare section SEC1 or SEC2, corresponding to {y=0, v_y>0}
  or {y=0, v_y<0}.
  Let $x$ be a point in the section, and suppose that flow at $x$ is
  transversal to the section.
  Compute the n-th iterate of the inverse Poincare map $P^{-n}(x)$ of the RTBP.
  This procedure also computes, as a side product, the integration time to
  intersect the Poincare section "n" times.

  \param[in] mu mass parameter for the RTBP
  \param[in] sec type of Poincare section (sec = SEC1 or SEC2).

  \param[in] cuts 
  number of iterates of the Poincare map (cuts=n: n cuts with the Poincare
  section).

  \param[in,out] x 
  Initial point, 4 coordinates: (X, Y, P_X, P_Y). 
  On return of the this function, it holds the image point $P^{-n}(x)$.

  \param[out] ti
  On return, it holds the integration time to intersect the Poincare section
  "n" times. Since we are computing the inverse Poincare map, this
  integration time "ti" must be negative.

  \return
  Returns a non-zero error code to indicate an error and 0 to indicate
  success.
  If an integration error is encountered, the function returns a non-zero
  value.

  \remark
  We are only interested in orbits that go around the origin. For some values
  of the energy, the orbit has "contractible" loops, which we want to avoid.
  For this, we impose that TWO CONSECUTIVE ITERATES do NOT lie both to the
  right or to the left of the origin.

  \remark
  On successful return of this function, the point $x$ is exactly on the
  section, i.e. we set coordinate $y$ exactly equal to zero.
*/    
int prtbp_inv(double mu, section_t sec, int cuts, double x[DIM], double *ti);

#endif // PRTBP_H_INCLUDED
