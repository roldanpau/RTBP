// ====================================================================
// Derivative of Poincare map of the Restricted Three Body Problem (2D)
// ====================================================================
// FILE:          $RCSfile: dprtbp_2d.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2013-01-25 12:07:06 $
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
#include <frtbp.h>	// DIMV
#include <prtbp.h>	// prtbp
#include <hinv.h>
#include <math.h>	// fabs
#include "dprtbp.h"	// dprtbp
#include "dprtbp_2d.h"	// ERR_IFT

const int ERR_VECTORFIELD=2;
const int ERR_DPRTBP=3;
const int ERR_POINCARE=4;

// name OF FUNCTION: dprtbp_2d
// CREDIT: 
//
// PURPOSE
// =======
// Consider the RTBP. Let "sec" be the Poincare section.
// Fix the Hamiltonian to a given value "H". 
// Using this energy condition, we can work with only two variables, $(x,p_x)$.
// The third variable $y$ is determined by the Poincare section, and the
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
// sec
//    type of Poincare section = {SEC1,SEC2}
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
// CALLS TO: hinv, prtbp, rtbp, dprtbp, set_dprtbp_2d

int dprtbp_2d(double mu, section_t sec, double H, int cuts, double p[2], double dp2d[4])
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
   hinv(mu,sec,H,x);

   // Compute derivative of Poincare map.
   if(dprtbp(mu,sec,cuts,x,dp))
   {
      fprintf(stderr, \
	    "dprtbp_2d: error computing derivative of poincare map\n");
      return(ERR_DPRTBP);
   }

   // Vectorfield evaluated at point "x"
   if(rtbp(0.0,x,f,&mu))
   {
      fprintf(stderr, "dprtbp_2d: error computing vectorfield\n");
      return(ERR_VECTORFIELD);
   }

   // Vectorfield evaluated at point "y=P(x)"
   y[0]=x[0]; y[1]=x[1]; y[2]=x[2]; y[3]=x[3];
   if(prtbp(mu,sec,cuts,y,&ti))
   {
      fprintf(stderr, "dprtbp_2d: error computing poincare map\n");
      return(ERR_POINCARE);
   }
   if(rtbp(0.0,y,g,&mu))
   {
      fprintf(stderr, "dprtbp_2d: error computing vectorfield");
      exit(ERR_VECTORFIELD);
   }

   // Set the 2d derivative $DP(p)$.
   return(set_dprtbp_2d(dp,f,g,dp2d));
}

// name OF FUNCTION: dprtbp_2d_inv
// CREDIT: 
//
// PURPOSE
// =======
// Consider the RTBP. Let "sec" be the Poincare section.
// Fix the Hamiltonian to a given value "H". 
// Using this energy condition, we can work with only two variables, $(x,p_x)$.
// The third variable $y$ is determined by the Poincare section, and the
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
// sec
//    type of Poincare section = {SEC1,SEC2}
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
// CALLS TO: hinv, prtbp, rtbp, dprtbp, set_dprtbp_2d

int dprtbp_2d_inv(double mu, section_t sec, double H, int cuts, double p[2], double dp2d[4])
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
   hinv(mu,sec,H,x);

   // Compute derivative of inverse Poincare map.
   if(dprtbp_inv(mu,sec,cuts,x,dp))
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
   if(prtbp_inv(mu,sec,cuts,y,&ti))
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

// name OF FUNCTION: set_dprtbp_2d
// CREDIT: 
//
// PURPOSE
// =======
// Compute the 2D derivative $DP(x,p_x)$ from the 4D variationals
// $DP(x,y,p_x,p_y)$, and from the vectorfield of the RTBP evaluated at the
// point $p=(x,0,p_x,p_y)$ and $\phi_T(p)$.
//
// Let 
// \[ dp=(u11 & u12 & u13 & u14\\ 
//        u21 & u22 & u23 & u24\\ 
//        u31 & u32 & u33 & u34\\ 
//        u41 & u42 & u43 & u44)
// \]
// be the 4D variationals.
// Let $f$ be the vectorfield at $p=(x,0,p_x,p_y)$.
// Let $g$ be the vectorfield at $\phi_T(p)$.
// 
// Then, the 2D derivative $dp2d$ is
// \[ dp2d=(D11 & D12\\
//          D21 & D22),
// \]
// where
// \[ D11 = u11 + u14*f3/f2 + g1*(-1.0/g2*(u21+u24*f3/f2)), 
//    D12 = u13 - u14*f1/f2 + g1*(-1.0/g2*(u23-u24*f1/f2)), 
//    D21 = u31 + u34*f3/f2 + g3*(-1.0/g2*(u21+u24*f3/f2)), 
//    D22 = u33 - u34*f1/f2 + g3*(-1.0/g2*(u23-u24*f1/f2)).
// \]
// See the notes ../doc/splitting.pdf for the derivation of these formulae.
//
// PARAMETERS
// ==========
// dp
//    4D variationals.
// 
// f
//    Vectorfield at $p=(x,0,p_x,p_y)$.
//
// g
//    Vectorfield at $\phi_T(p)$.
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
// ERR_IFT
//    Cannot apply Implicit Function Theorem to function $p_y(x,p_x)$.
//
// NOTES
// =====
// In order to apply the Implicit Function Theorem to express the function
// $p_y=p_y(x,p_x)$ as a graph, we need that
// \[ \partial H/\partial p_y = \dot y \neq 0, \]
// that is, the second component of the vector field must be nonzero.

const int ERR_IFT=1;

int set_dprtbp_2d(double dp[DIMV], double f[DIM], double g[DIM], double
      dp2d[4])
{
   double f1,f2,f3;
   double g1,g2,g3;
   double u11, u13, u14, u21, u23, u24, u31, u33, u34;

   // Auxiliary variables
   double dpydx, dTdx;
   double dpydpx, dTdpx;

   // Check that the second component of the vector field is nonzero.
   if(fabs(f[1])<1.e-10)
   {
      fprintf(stderr, \
	    "dprtbp_2d: can't apply IFT to function py(x,px), " \
	    "because y_dot = %le is too small\n", f[1]);
      return(ERR_IFT);
   }
   f1=f[0]; f2=f[1]; f3=f[2];
   g1=g[0]; g2=g[1]; g3=g[2];
   u11=dp[0]; u13=dp[2]; u14=dp[3];
   u21=dp[4]; u23=dp[6]; u24=dp[7];
   u31=dp[8]; u33=dp[10]; u34=dp[11];
   dpydx=f3/f2; 
   dTdx=-1.0/g2*(u21+u24*dpydx);
   dpydpx=-f1/f2; 
   dTdpx=-1.0/g2*(u23+u24*dpydpx);
   dp2d[0] = u11+u14*dpydx+g1*dTdx;
   dp2d[1] = u13+u14*dpydpx+g1*dTdpx;
   dp2d[2] = u31+u34*dpydx+g3*dTdx;
   dp2d[3] = u33+u34*dpydpx+g3*dTdpx;
   return(0);
}
