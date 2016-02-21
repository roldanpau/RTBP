// ======================================================
// Poincare map of the Restricted Three Body Problem (2D)
// ======================================================
// FILE:          $RCSfile: sec1sec2_module.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2012-12-20 11:02:51 $
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
// sec1sec2
//    This function takes a point $p$ on section S1 to section S2 by forward
//    integration.
// sec1sec2_inv
//    This function takes a point $p$ on section S1 to section S2 by backward
//    integration.

#include <stdio.h>	// fprintf
#include <stdlib.h>	// EXIT_FAILURE
#include <math.h>	// sqrt, fabs
#include <rtbp.h>	// DIM, rtbp
#include <hinv.h>
#include <prtbp.h>

// name OF FUNCTION: sec1sec2
//
// PURPOSE
// =======
// Consider the RTBP. 
// Let S1 be the Poincare section {y=0, vy>0}.
// Let S2 be the Poincare section {y=0, vy<0}.
// Fix the Hamiltonian to a given value "H". 
// Using this energy condition, we can work with only two variables, $(x,p_x)$.
// The third variable $y$ is 0 since we look at the Poincare section, and the
// fourth variable $p_y$ can be obtained from the energy condition.
// This function takes a point $p$ on section S1 to section S2 by forward
// integration.
//
// PARAMETERS
// ==========
// mu
//    mass parameter for the RTBP
// H
//    energy value
// p
//    Initial point, 2 coordinates: p=(x,p_x). 
//    On return of the this function, it holds the image point on the section
//    S2.
// ti
//    On exit, *ti holds the integration time to reach the section S2.
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
// CALLED FROM:

int sec1sec2(double mu, double H, double p[2], double *ti)
{
   double x[DIM];
   x[0]=p[0];	// x
   x[1]=0.0;	// y
   x[2]=p[1]; 	// p_x

   // Compute x[3]=p_y by inverting the Hamiltonian.
   // Recall that p is on section SEC1.
   if(hinv(mu,SEC1,H,x))
   {
      fprintf(stderr, "sec1sec2: error inverting the Hamiltonian\n");
      return(1);
   }

   // Compute first intersection with section S2 (and integration time).
   // On exit, "prtbp" guarantees that x is exactly on Poincare section S2.
   if(prtbp(mu,SEC2,1,x,ti))
   {
      fprintf(stderr, "sec1sec2: error computing poincare map\n");
      return(1);
   }
   // Set the image point on section S2.
   p[0]=x[0]; p[1]=x[2];
   return(0);
}

int sec1sec2_inv(double mu, double H, double p[2], double *ti)
{
   double x[DIM];
   x[0]=p[0];	// x
   x[1]=0.0;	// y
   x[2]=p[1]; 	// p_x

   // Compute x[3]=p_y by inverting the Hamiltonian.
   // Recall that p is on section SEC1.
   if(hinv(mu,SEC1,H,x))
   {
      fprintf(stderr, "sec1sec2_inv: error inverting the Hamiltonian\n");
      return(1);
   }

   // Compute first intersection with section S2 (and integration time).
   // On exit, "prtbp_inv" guarantees that x is exactly on Poincare section S2.
   if(prtbp_inv(mu,SEC2,1,x,ti))
   {
      fprintf(stderr, "sec1sec2_inv: error computing inverse poincare map\n");
      return(1);
   }
   // Set the image point on section S2.
   p[0]=x[0]; p[1]=x[2];
   return(0);
}
