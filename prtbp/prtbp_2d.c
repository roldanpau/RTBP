// ======================================================
// Poincare map of the Restricted Three Body Problem (2D)
// ======================================================
// FILE:          $RCSfile: prtbp_2d.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2012-12-19 10:58:59 $
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
#include <hinv.h>	// hinv
#include "prtbp.h"	// section_t, prtbp, prtbp_inv

// name OF FUNCTION: prtbp_2d
// CREDIT: 
//
// PURPOSE
// =======
// Consider the RTBP. Let "sec" be the Poincare section.
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
// sec
//    type of Poincare section (sec = SEC1 or SEC2).
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
// CALLS TO: hinv, prtbp
//
// CALLED FROM: main

int prtbp_2d(double mu, section_t sec, double H, int cuts, double p[2], double *ti)
{
   double x[DIM];

   // auxiliary variables
   int result;

   x[0]=p[0];	// x
   x[1]=0.0;	// y
   x[2]=p[1]; 	// p_x

   // Compute x[3]=p_y by inverting the Hamiltonian.
   result = hinv(mu,sec,H,x);
   if(result)
   {
      fprintf(stderr, "prtbp_2d: error inverting the Hamiltonian\n");
      return(1);
   }

   // Compute Poincare map (and integration time).
   // On exit, "prtbp" guarantees that x is precisely on the Poincare section.
   if(prtbp(mu,sec,cuts,x,ti))
   {
      fprintf(stderr, "prtbp_2d: error computing poincare map\n");
      return(1);
   }
   // Set the image point $P^n(p)$.
   p[0]=x[0]; 	// x'
   p[1]=x[2];	// p_x'
   return(0);
}

// name OF FUNCTION: prtbp_2d_inv
// CREDIT: 
//
// PURPOSE
// =======
// Consider the RTBP. Let "sec" be the Poincare section.
// Fix the Hamiltonian to a given value "H". 
// Using this energy condition, we can work with only two variables, $(x,p_x)$.
// The third variable $y$ is 0 since we look at the Poincare section, and the
// fourth variable $p_y$ can be obtained from the energy condition.
// This function computes the n-th iterate $P^{-n}$ of the inverse 2D Poincare map:
//    (x',px') = P^{-n}(x,px).
//
// PARAMETERS
// ==========
// mu
//    mass parameter for the RTBP
// sec
//    type of Poincare section (sec = SEC1 or SEC2).
// H
//    energy value
// cuts
//    number of iterates of the Poincare map (cuts=n: n cuts with the
//    Poincare section).
// p
//    Initial point, 2 coordinates: p=(x,p_x). 
//    On return of the this function, it holds the image point $P^{-n}(p)$.
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
// CALLS TO: hinv, prtbp_inv
//
// CALLED FROM: main

int prtbp_2d_inv(double mu, section_t sec, double H, int cuts, double p[2], double *ti)
{
   double x[DIM];

   // auxiliary variables
   int result;

   x[0]=p[0];	// x
   x[1]=0.0;	// y
   x[2]=p[1]; 	// p_x

   // Compute x[3]=p_y by inverting the Hamiltonian.
   result = hinv(mu,sec,H,x);
   if(result)
   {
      fprintf(stderr, "prtbp_2d_inv: error inverting the Hamiltonian\n");
      return(1);
   }

   // Compute inverse Poincare map (and integration time).
   // On exit, "prtbp_inv" guarantees that x is precisely on the Poincare section.
   if(prtbp_inv(mu,sec,cuts,x,ti))
   {
      fprintf(stderr, "prtbp_2d_inv: error computing inverse poincare map\n");
      return(1);
   }
   // Set the image point $P^{-n}(p)$.
   p[0]=x[0]; 	// x'
   p[1]=x[2];	// p_x'
   return(0);
}
