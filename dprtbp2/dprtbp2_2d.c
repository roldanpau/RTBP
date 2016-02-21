// ====================================================================
// Derivative of Poincare map of the Restricted Three Body Problem (2D)
// ====================================================================
// FILE:          $RCSfile: dprtbp2_2d.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2011-03-18 14:08:30 $
//
// FUNCTIONS
// =========
//
// dprtbp_2d
// ---------
// Compute the derivative $DP$ of the 2D Poincare map.
//
// dprtbp_2d_inv
// -------------
// Compute the derivative $DP^{-1}$ of the inverse 2D Poincare map.
//
// set_dprtbp_2d
// -------------
// Compute the 2D derivative $DP(x,p_x)$ from the 4D variationals
// $DP(x,y,p_x,p_y)$, and from the vectorfield of the RTBP evaluated at the
// point $p=(x,0,p_x,p_y)$ and $\phi_T(p)$.

#include <stdio.h>	// fprintf
#include <stdlib.h>	// EXIT_FAILURE
#include <rtbp.h>	// DIM, rtbp
#include <prtbp2.h>
#include <hinv2.h>
#include <math.h>	// fabs
#include "dprtbp2.h"
#include "dprtbp_2d.h"	// ERR_IFT, ERR_VECTORFIELD, ..., set_dprtbp_2d

// name OF FUNCTION: dprtbp2_2d
// CREDIT: 
//
// PURPOSE
// =======
// Consider the RTBP. Let S be the Poincare section {y=0}.
// Fix the Hamiltonian to a given value "H". 
// Using this energy condition, we can work with only two variables, $(x,p_x)$.
// The third variable $y$ is 0 since we look at the Poincare section, and the
// fourth variable $p_y$ can be obtained from the energy condition.
// This function computes the derivative $DP(x,p_x)$ of the 2D Poincare map.
// The 2D derivative is obtained from the variationals (i.e. the full 4D
// derivative) in the auxiliary function "set_dprtbp_2d".
//
// PARAMETERS
// ==========
// mu
//    mass parameter for the RTBP
//
// H
//    energy value
//
// cuts
//    number of iterates of the Poincare map (cuts=n: n cuts with the
//    Poincare section).
// 
// p
//    Argument of the 2d derivative, 2 coordinates: p=(x,p_x). 
//
// dp2d
//    On exit, dp2d holds the (2-by-2) derivative $DP$ of the 2D Poincare
//    map.
// 
// RETURN VALUE
// ============
// Returns a non-zero error code to indicate an error and 0 to indicate
// success.
// If an error is encountered, the function returns a non-zero value:
//
// ERR_VECTORFIELD
//    Error computing RTBP vectorfield (collision).
//
// ERR_DPRTBP
//    Error computing derivative of 4D Poincare map of RTBP.
//
// ERR_POINCARE
//    Error computing Poincare map $P(x,p_x)$. If we can't compute the
//    function $P$, we can't compute its derivative!
//
// ERR_IFT
//    Cannot apply Implicit Function Theorem to function $p_y(x,p_x)$.
//
// NOTES
// =====
// In order to apply the Implicit Function Theorem to express the function
// $p_y=p_y(x,p_x)$ as a graph, we need that
// \[ \partial H/\partial p_y = \dot y \neq 0, \]
// that is, the second component of the vector field must be nonzero.
//
// CALLS TO: hinv2, prtbp2, rtbp, dprtbp2, set_dprtbp_2d

int dprtbp2_2d(double mu, double H, int cuts, double p[2], double dp2d[4])
{
   double x[DIM];
   double y[DIM];	// y = P(x) image of "x" under Poincare map
   double f[DIM];	// 4D vectorfield evaluated at "x"
   double g[DIM];	// 4D vectorfield evaluated at "y"

   // 4D derivative of poincare map
   double dp[DIMV];

   // Auxiliary variables
   double ti;

   x[0]=p[0]; x[1]=0.0;	// x, y
   x[2]=p[1]; 		// p_x

   // Compute x[3]=p_y by inverting the Hamiltonian.
   hinv2(mu,H,x);

   // Compute derivative of Poincare map.
   if(dprtbp2(mu,cuts,x,dp))
   {
      fprintf(stderr, \
	    "dprtbp2_2d: error computing derivative of poincare map\n");
      return(ERR_DPRTBP);
   }

   // Vectorfield evaluated at point "x"
   if(rtbp(0.0,x,f,&mu))
   {
      fprintf(stderr, "dprtbp2_2d: error computing vectorfield\n");
      return(ERR_VECTORFIELD);
   }

   // Vectorfield evaluated at point "y=P(x)"
   y[0]=x[0]; y[1]=x[1]; y[2]=x[2]; y[3]=x[3];
   if(prtbp2(mu,cuts,y,&ti))
   {
      fprintf(stderr, "dprtbp2_2d: error computing poincare map\n");
      return(ERR_POINCARE);
   }
   if(rtbp(0.0,y,g,&mu))
   {
      fprintf(stderr, "dprtbp2_2d: error computing vectorfield");
      exit(ERR_VECTORFIELD);
   }

   // Set the 2d derivative $DP(p)$.
   return(set_dprtbp_2d(dp,f,g,dp2d));
}

// name OF FUNCTION: dprtbp2_2d_inv
// CREDIT: 
//
// PURPOSE
// =======
// Consider the RTBP. Let S be the Poincare section {y=0}.
// Fix the Hamiltonian to a given value "H". 
// Using this energy condition, we can work with only two variables, $(x,p_x)$.
// The third variable $y$ is 0 since we look at the Poincare section, and the
// fourth variable $p_y$ can be obtained from the energy condition.
// This function computes the derivative $DP(x,p_x)$ of the 2D Poincare map.
// The 2D derivative is obtained from the variationals (i.e. the full 4D
// derivative) in the auxiliary function "set_dprtbp_2d".
//
// PARAMETERS
// ==========
// mu
//    mass parameter for the RTBP
//
// H
//    energy value
//
// cuts
//    number of iterates of the Poincare map (cuts=n: n cuts with the
//    Poincare section).
// 
// p
//    Argument of the 2d derivative, 2 coordinates: p=(x,p_x). 
//
// dp2d
//    On exit, dp2d holds the (2-by-2) derivative $DP$ of the 2D Poincare
//    map.
// 
// RETURN VALUE
// ============
// Returns a non-zero error code to indicate an error and 0 to indicate
// success.
// If an error is encountered, the function returns a non-zero value:
//
// ERR_VECTORFIELD
//    Error computing RTBP vectorfield (collision).
//
// ERR_DPRTBP
//    Error computing derivative of 4D Poincare map of RTBP.
//
// ERR_POINCARE
//    Error computing Poincare map $P(x,p_x)$. If we can't compute the
//    function $P$, we can't compute its derivative!
//
// ERR_IFT
//    Cannot apply Implicit Function Theorem to function $p_y(x,p_x)$.
//
// NOTES
// =====
// In order to apply the Implicit Function Theorem to express the function
// $p_y=p_y(x,p_x)$ as a graph, we need that
// \[ \partial H/\partial p_y = \dot y \neq 0, \]
// that is, the second component of the vector field must be nonzero.
//
// CALLS TO: hinv2, prtbp2, rtbp, dprtbp2, set_dprtbp_2d

int dprtbp2_2d_inv(double mu, double H, int cuts, double p[2], double dp2d[4])
{
   double x[DIM];
   double y[DIM];	// y = P^{-1}(x) image of "x" under inverse Poincare map
   double f[DIM];	// 4D vectorfield evaluated at "x"
   double g[DIM];	// 4D vectorfield evaluated at "y"

   // 4D derivative of inverse poincare map
   double dp[DIMV];

   // Auxiliary variables
   double ti;

   x[0]=p[0]; x[1]=0.0;	// x, y
   x[2]=p[1];		// p_x

   // Compute x[3]=p_y by inverting the Hamiltonian.
   hinv2(mu,H,x);

   // Compute derivative of inverse Poincare map.
   if(dprtbp2_inv(mu,cuts,x,dp))
   {
      fprintf(stderr, \
	    "dprtbp_2d_inv: error computing derivative of poincare map\n");
      return(ERR_DPRTBP);
   }

   // Negative vectorfield evaluated at point "x"
   if(rtbp_inv(0.0,x,f,&mu))
   {
      fprintf(stderr, "dprtbp_2d_inv: error computing vectorfield\n");
      return(ERR_VECTORFIELD);
   }

   // Negative vectorfield evaluated at point "y=P^{-1}(x)"
   y[0]=x[0]; y[1]=x[1]; y[2]=x[2]; y[3]=x[3];
   if(prtbp2_inv(mu,cuts,y,&ti))
   {
      fprintf(stderr, "dprtbp_2d_inv: error computing poincare map\n");
      return(ERR_POINCARE);
   }
   if(rtbp_inv(0.0,y,g,&mu))
   {
      fprintf(stderr, "dprtbp_2d_inv: error computing vectorfield");
      exit(ERR_VECTORFIELD);
   }

   // Set the 2d derivative $DP(p)$.
   return(set_dprtbp_2d(dp,f,g,dp2d));
}

