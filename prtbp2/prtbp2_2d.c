// ======================================================
// Poincare map of the Restricted Three Body Problem (2D)
// ======================================================
// FILE:          $RCSfile: prtbp2_2d.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2012-12-20 11:23:15 $
//
// PURPOSE:
//
// NOTES:
//
// OVERALL METHOD:
//
// FUNCTIONS
// =========
//
// prtbp_2d
//    Compute the n-th iterate $P^n$ of the 2D Poincare map:
//       (x',px') = P^n(x,px)
//    in the RTBP.
// prtbp_2d_inv
//    Compute the inverse $P^{-1}$ of the 2D Poincare map:
//       (x',px') = P^{-1}(x,px)
//    in the RTBP.

#include <stdio.h>	// fprintf
#include <stdlib.h>	// EXIT_FAILURE
#include <math.h>	// sqrt, fabs
#include <rtbp.h>	// DIM, rtbp
#include <hinv2.h>
#include "prtbp2.h"

// name OF FUNCTION: prtbp2_2d
// CREDIT: 
//
// PURPOSE
// =======
// Consider the RTBP. Let S be the Poincare section {y=0, px<0}.
// Fix the Hamiltonian to a given value "H". 
// Using this energy condition, we can work with only two variables, $(x,p_x)$.
// The third variable $y$ is 0 since we look at the Poincare section, and the
// fourth variable $p_y$ can be obtained from the energy condition.
// This function computes the n-th iterate $P^n$ of the 2D Poincare map:
//    (x',px') = P^n(x,px).
//
// PARAMETERS
// ==========
// mu
//    mass parameter for the RTBP
// H
//    energy value
// cuts
//    number of iterates of the Poincare map (cuts=n: n cuts with the
//    Poincare section).
// p
//    Initial point, 2 coordinates: p=(x,p_x). 
//    On return of the this function, it holds the image point $P^n(p)$.
// ti
//    On exit, *ti holds the integration time to reach the Poincare section.
// 
// RETURN VALUE
// ============
// Returns a non-zero error code to indicate an error and 0 to indicate
// success.
//
// NOTES
// =====
//
// CALLS TO: hinv2, prtbp2
//
// CALLED FROM: main

int prtbp2_2d(double mu, double H, int cuts, double p[2], double *ti)
{
   double x[DIM];
   x[0]=p[0];	// x
   x[1]=0.0;	// y
   x[2]=p[1]; 	// p_x

   // Compute x[3]=p_y by inverting the Hamiltonian.
   if(hinv2(mu,H,x))
   {
      fprintf(stderr, "prtbp2_2d: error inverting the Hamiltonian\n");
      return(1);
   }

   // Compute Poincare map (and integration time).
   // On exit, "prtbp" guarantees that x is exactly on Poincare section.
   if(prtbp2(mu,cuts,x,ti))
   {
      fprintf(stderr, "prtbp2_2d: error computing poincare map\n");
      return(1);
   }
   // Set the image point $P^n(p)$.
   p[0]=x[0]; p[1]=x[2];
   return(0);
}

// name OF FUNCTION: prtbp_2d_inv
// CREDIT: 
//
// PURPOSE
// =======
// Consider the RTBP. Let S be the Poincare section {y=0}.
// Fix the Hamiltonian to a given value "H". 
// Using this energy condition, we can work with only two variables, $(x,p_x)$.
// The third variable $y$ is 0 since we look at the Poincare section, and the
// fourth variable $p_y$ can be obtained from the energy condition.
// This function computes the inverse $P^{-1}$ of the 2D Poincare map:
//    (x',px') = P^{-1}(x,px).
//
// PARAMETERS
// ==========
// mu
//    mass parameter for the RTBP
// H
//    energy value
// cuts
//    number of iterates of the Poincare map (cuts=n: n cuts with the
//    Poincare section).
// p
//    Initial point, 2 coordinates: p=(x,p_x). 
//    On return of the this function, it holds the image point $P^{-1}(p)$.
// ti
//    On exit, *ti holds the integration time to reach the Poincare section.
//    Since we are computing the inverse Poincare map, the integration time
//    "ti" must be negative.
// 
// RETURN VALUE
// ============
// Returns a non-zero error code to indicate an error and 0 to indicate
// success.
//
// NOTES
// =====
//
// CALLS TO: hinv2, prtbp2_inv
//
// CALLED FROM: main

int prtbp2_2d_inv(double mu, double H, int cuts, double p[2], double *ti)
{
   double x[DIM];
   x[0]=p[0];	// x
   x[1]=0.0;	// y
   x[2]=p[1]; 	// p_x

   // Compute x[3]=p_y by inverting the Hamiltonian.
   if(hinv2(mu,H,x))
   {
      fprintf(stderr, "prtbp_2d_inv: error inverting the Hamiltonian\n");
      return(1);
   }

   // Compute inverse Poincare map (and integration time).
   // On exit, "prtbp_inv" guarantees that x is exactly on Poincare section.
   if(prtbp2_inv(mu,cuts,x,ti))
   {
      fprintf(stderr, "prtbp_2d_inv: error computing inverse poincare map\n");
      return(1);
   }
   // Set the image point $P^{-1}(p)$.
   p[0]=x[0]; p[1]=x[2];
   return(0);
}
