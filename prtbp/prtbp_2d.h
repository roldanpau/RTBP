/*! \file
    \brief Poincare map of the Restricted Three Body Problem (2D)

    $Author: roldan $
    $Date: 2013-03-26 22:22:14 $
*/

#ifndef PRTBP_2D_H_INCLUDED
#define PRTBP_2D_H_INCLUDED

#include "prtbp.h"	// section_t

/**
  Compute the n-th iterate $P^n$ of the 2D Poincare map:
  \f[
     (x',px') = P^n(x,px)
  \f]
  in the RTBP.

  Consider the RTBP in rotating coordinates.
  Let S be the Poincare section SEC1 or SEC2, corresponding to {y=0, v_y>0}
  or {y=0, v_y<0}.
 Fix the Hamiltonian to a given value "H". 
 Using this energy condition, we can work with only two variables, $(x,p_x)$.
 The third variable $y$ is 0 since we look at the Poincare section, and the
 fourth variable $p_y$ can be obtained from the energy condition.
 This function computes the n-th iterate $P^n$ of the 2D Poincare map:
 \f[
    (x',px') = P^n(x,px).
 \f]
  This procedure also computes, as a side product, the integration time to
  intersect the Poincare section "n" times.

  \param[in] mu mass parameter for the RTBP
  \param[in] sec type of Poincare section (sec = SEC1 or SEC2).

  \param[in] H energy value
  \param[in] cuts 
  number of iterates of the Poincare map (cuts=n: n cuts with the Poincare
  section).

  \param[in,out] p 
  Initial point, 2 coordinates: p=(x,p_x). 
  On return of the this function, it holds the image point $P^n(p)$.

  \param[out] ti
  On return, it holds the integration time to intersect the Poincare
  section "n" times.

  \return
  Returns a non-zero error code to indicate an error and 0 to indicate
  success.
*/    

int prtbp_2d(double mu, section_t sec, double H, int cuts, double p[2], double *ti);

/**
  Compute the n-th iterate $P^{-n}$ of the inverse 2D Poincare map:
  \f[
     (x',px') = P^{-n}(x,px)
  \f]
  in the RTBP.

  Consider the RTBP in rotating coordinates.
  Let S be the Poincare section SEC1 or SEC2, corresponding to {y=0, v_y>0}
  or {y=0, v_y<0}.
 Fix the Hamiltonian to a given value "H". 
 Using this energy condition, we can work with only two variables, $(x,p_x)$.
 The third variable $y$ is 0 since we look at the Poincare section, and the
 fourth variable $p_y$ can be obtained from the energy condition.
 This function computes the n-th iterate $P^{-n}$ of the inverse 2D Poincare map:
 \f[
    (x',px') = P^{-n}(x,px).
 \f]
  This procedure also computes, as a side product, the integration time to
  intersect the Poincare section "n" times.

  \param[in] mu mass parameter for the RTBP
  \param[in] sec type of Poincare section (sec = SEC1 or SEC2).

  \param[in] H energy value
  \param[in] cuts 
  number of iterates of the Poincare map (cuts=n: n cuts with the Poincare
  section).

  \param[in,out] p 
  Initial point, 2 coordinates: p=(x,p_x). 
  On return of the this function, it holds the image point $P^{-n}(p)$.

  \param[out] ti
  On return, it holds the integration time to intersect the Poincare
  section "n" times.
  Since we are computing the inverse Poincare map, the integration time
  "ti" must be negative.

  \return
  Returns a non-zero error code to indicate an error and 0 to indicate
  success.
*/    

int prtbp_2d_inv(double mu, section_t sec, double H, int cuts, double p[2], double *ti);

#endif // PRTBP_2D_H_INCLUDED
