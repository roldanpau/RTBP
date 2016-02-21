/*! \file
    \brief Flow of the RTBP in Delaunay coords

    $Author: roldan $
    $Date: 2013-03-26 22:18:08 $
*/

#include "rtbp.h"	// DIM

/// dimension of the (planar) RTBP variational equations
#define DIMV 16

/**
  Flow of the RTBP in Delaunay coords.

  Compute the flow $\phi(t,x)$ of the RTBP for a given time and initial
  condition. 
  The time $t$ may be positive or negative, allowing for forward or backward
  integration.
  The trajectory is integrated numerically.
 
  \param[in] mu
     mass parameter for the RTBP

  \param[in] t1
     integration time 

  \param[in,out] x[DIM]
     Initial condition, 4 coordinates: (l, L, g, G). 
     On return of this function, it holds the final point $\phi(t,x)$.
  
  \return
  a non-zero error code to indicate an error and 0 to indicate success.
  If an integration error is encountered, the function returns a non-zero
  value.
 
  \remark
  A Runge-Kutta Prince-Dormand (8,9) method is used to solve the ODE. 
  We request absolute local error 10^{-16} and relative local error 0.
  
 */

int frtbp_del(double mu, double t1, double x[DIM]);
