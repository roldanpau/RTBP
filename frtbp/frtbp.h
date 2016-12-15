/*! \file
    \brief Flow of the Restricted Three Body Problem

    $Author: roldan $
    $Date: 2013-03-26 22:17:52 $
*/

/*! \mainpage Restricted Three Body Problem Library

  \section Introduction

  This is a library of functions for different computations related to the
  Restricted Three Body Problem, or RTBP for short. The software is available
  upon request (just email me), so that I can keep track of who's using it
  and why. 

  This software has been developed for my own research, in particular for the
  following paper:
  "Diffusion along mean motion resonance in the restricted planar three-body
  problem" by J.Fejoz, M.Guardia, V.Kaloshin and P.Roldan.
  (<a href="http://arxiv.org/abs/1109.2892">link</a>)

  The kind of functions that you will find here are
  - Numerical integration of the RTBP equations (in Cartesian variables)
  - Numerical integration of the RTBP equations (in Delaunay variables)
  - Cartesian to Delaunay change of coordinates
  - Numerical integration of the RTBP variational equations 
  - Poincare section and associated Poincare map (in Cartesian)
  - Poincare section and associated Poincare map (in Delaunay)
  - Numerical computation of (some family of) periodic orbits
  - Hyperbolic splitting (eigenvalues/vects) associated to a fixed point of
  the Poincare map
  - Numerical computation of hyperbolic invariant manifolds associated to
  fixed point
  - Homoclinic intersection of invariant manifolds
  - Numerical computation of splitting angle of the manifolds at the
  homoclinic intersection

  \section Acknowledgements

  A big thanks goes to Angel Jorba for altruistically sharing his ideas about
  RTBP computations with me. Moreover, he provided me with the high-precision 
  Taylor method integrator which is at the base of some of my computations.
  See 
  <a href="http://www.maia.ub.edu/~angel/soft.html">http://www.maia.ub.edu/~angel/soft.html</a>.
 */

#include "rtbp.h"

/// dimension of the (planar) RTBP variational equations
#define DIMV 16

/**
  Derivative of the flow of the RTBP.

  Compute the derivative of the flow \f$D\phi(t,x)\f$ with respect to $x$ of the
  RTBP for a given time and initial condition. 
  The time $t$ may be positive or negative, allowing for forward or backward
  integration.
  This is done by numerically integrating the variational equations.
 
  \param[in] mu_loc
     mass parameter for the RTBP

  \param[in] t1
     integration time 

  \param[in] x
     Argument of the derivative, 4 coordinates: \f$ (X, Y, P_X, P_Y). \f$

  \param[out] dphi
     On return of this function, it holds the derivative \f$D\phi(t,x)\f$ with
     respect to $x$.
  
  \return
  a non-zero error code to indicate an error and 0 to indicate success.
 
  \remark
  Backwards integration of the variational equations is used to obtain the
  derivative of the inverse Poincare map (see \ref dprtbp).

  \remark
  A Taylor method (provided by Angel Jorba) is used to solve the variational
  equations.
 
  \remark
  We follow the convention to place the large mass (Sun) to the LEFT of the
  origin, and the small mass (Jupiter) to the RIGHT.
  
 */

int dfrtbp(double mu_loc, double t1, double x[DIM], double dphi[DIMV]);

/**
  Flow of the Restricted Three Body Problem

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
     Initial condition, 4 coordinates: \f$ (X, Y, P_X, P_Y). \f$
     On return of this function, it holds the final point $\phi(t,x)$.
  
  \return
  a non-zero error code to indicate an error and 0 to indicate success.
  If an integration error is encountered, the function returns a non-zero
  value.
 
  \remark
  A Taylor method (provided by Angel Jorba) is used to solve the ODE.
 
  \remark
  We follow the convention to place the large mass (Sun) to the LEFT of the
  origin, and the small mass (Jupiter) to the RIGHT.
  
 */

int frtbp(double mu, double t1, double x[DIM]);
