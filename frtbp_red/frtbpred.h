/*! \file
    \brief Flow of the Reduced Restricted Three Body Problem

    $Author: roldan $
    $Date: 2013-03-26 22:18:24 $
*/

#include "rtbpred.h"	// DIMRED

/**
  Flow of the Reduced Restricted Three Body Problem

  Compute the flow $\phi(s,x)$ of the reduced RTBP for a given time and
  initial condition. Reduced means that we identify $l$ with time (in
  frtbp_red_l) or $g$ with time (in frtbp_red_g) and then we do not need to
  integrate this variable.
  The time $s$ may be positive or negative, allowing for forward or backward
  integration.
  The trajectory is integrated numerically.
 
  \param[in] mu
     mass parameter for the RTBP

  \param[in] s1
     integration time 

  \param[in,out] x[DIMRED]
     Initial condition, 6 coordinates: (l, L, g, G, t, I). 
     On return of this function, it holds the final point $\phi(s,x)$.
  
  \return
  a non-zero error code to indicate an error and 0 to indicate success.
  If an integration error is encountered, the function returns a non-zero
  value.
 
  \remark
  A Runge-Kutta Prince-Dormand (8,9) method is used to solve the ODE. 
  We request absolute local error 10^{-16} and relative local error 0.
  */

int frtbp_red_l(double mu, double s1, double x[DIMRED]);
int frtbp_red_g(double mu, double s1, double x[DIMRED]);
