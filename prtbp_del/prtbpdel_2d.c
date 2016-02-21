// ======================================================
// Poincare map of the Restricted Three Body Problem (2D)
// ======================================================
// FILE:          $RCSfile: prtbpdel_2d.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2012-11-26 10:42:05 $
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
// prtbp_del_2d
//    Compute the n-th iterate $P^n$ of the 2D Poincare map:
//       (g',G') = P^n(g,G)
//    in the RTBP.
// prtbp_del_2d_inv
//    Compute the inverse $P^{-1}$ of the 2D Poincare map:
//       (g',G') = P^{-1}(g,G)
//    in the RTBP.
// sec1sec2_del
//    This function takes a point $p$ on section SEC1 to section SEC2 by forward
//    integration.
// sec1sec2_del_inv
//    This function takes a point $p$ on section SEC1 to section SEC2 by backward
//    integration.


#include <stdio.h>	// fprintf
#include <stdlib.h>	// EXIT_FAILURE

// Note: do not compile with flag -std=c99, or else M_PI will be undefined.
#include <math.h>	// M_PI, sqrt, fabs

#include <rtbp.h>	// DIM, rtbp
#include <hinvdel.h>	// hinv_del
#include "prtbpdel.h"	// section_t, prtbp_del, prtbp_del_inv

// name OF FUNCTION: prtbp_del_2d
// CREDIT: 
//
// PURPOSE
// =======
// Consider the RTBP. Let "sec" be the Poincare section.
// Fix the Hamiltonian to a given value "H". 
// Using this energy condition, we can work with only two variables, $(g,G)$.
// The third variable $l$ is 0 since we look at the Poincare section, and the
// fourth variable $L$ can be obtained from the energy condition.
// This function computes the n-th iterate $P^n$ of the 2D Poincare map:
//    (g',G') = P^n(g,G).
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
//    Initial point, 2 coordinates: p=(g,G). 
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
// CALLS TO: hinv_del, prtbp_del
//
// CALLED FROM: main

int prtbp_del_2d(double mu, section_t sec, double H, int cuts, double p[2], double *ti)
{
   double x[DIM];

   switch(sec)
   {
      case SEC1 :
	 {
	    x[0] = 0.0;		// l=0
	    break;
	 }
      case SEC2 :
	 {
	    x[0] = M_PI;	// l=pi
	    break;
	 }
   }
   x[2]=p[0];	// g
   x[3]=p[1]; 	// G

   // Compute x[1]=L by inverting the Hamiltonian.
   if(hinv_del(mu,H,x))
   {
      fprintf(stderr, "prtbp_del_2d: error inverting the Hamiltonian\n");
      return(1);
   }

   // Compute Poincare map (and integration time).
   // On exit, "prtbp_del" guarantees that x is extremely close to Poincare section.
   if(prtbp_del(mu,sec,cuts,x,ti))
   {
      fprintf(stderr, "prtbp_del_2d: error computing poincare map\n");
      return(1);
   }
   // Set the image point $P^n(p)$.
   p[0]=x[2]; 	// g
   p[1]=x[3];	// G
   return(0);
}

// name OF FUNCTION: prtbp_del_2d_inv
// CREDIT: 
//
// PURPOSE
// =======
// Consider the RTBP. Let "sec" be the Poincare section.
// Fix the Hamiltonian to a given value "H". 
// Using this energy condition, we can work with only two variables, $(g,G)$.
// The third variable $l$ is 0 since we look at the Poincare section, and the
// fourth variable $L$ can be obtained from the energy condition.
// This function computes the inverse $P^{-1}$ of the 2D Poincare map:
//    (g',G') = P^{-1}(g,G).
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
//    Initial point, 2 coordinates: p=(g,G). 
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
// CALLS TO: hinv_del, prtbp_del_inv
//
// CALLED FROM: main

int prtbp_del_2d_inv(double mu, section_t sec, double H, int cuts, double p[2], double *ti)
{
   double x[DIM];
   switch(sec)
   {
      case SEC1 :
	 {
	    x[0] = 0.0;		// l=0
	    break;
	 }
      case SEC2 :
	 {
	    x[0] = M_PI;	// l=pi
	    break;
	 }
   }
   x[2]=p[0];	// g
   x[3]=p[1]; 	// G

   // Compute x[1]=L by inverting the Hamiltonian.
   if(hinv_del(mu,H,x))
   {
      fprintf(stderr, "prtbp_del_2d_inv: error inverting the Hamiltonian\n");
      return(1);
   }

   // Compute inverse Poincare map (and integration time).
   // On exit, "prtbp_del_inv" guarantees that x is extremely close to Poincare section.
   if(prtbp_del_inv(mu,sec,cuts,x,ti))
   {
      fprintf(stderr, "prtbp_del_2d_inv: error computing inverse poincare map\n");
      return(1);
   }
   // Set the image point $P^{-1}(p)$.
   p[0]=x[2]; 	// g
   p[1]=x[3];	// G
   return(0);
}


// name OF FUNCTION: sec1sec2_del
//
// PURPOSE
// =======
// Consider the RTBP. 
// Let SEC1 be the Poincare section {l=0}.
// Let SEC2 be the Poincare section {l=pi}.
// Fix the Hamiltonian to a given value "H". 
// Using this energy condition, we can work with only two variables, $(g,G)$.
// The third variable is determined by the Poincare section, and the fourth
// variable $L$ can be obtained from the energy condition.
// This function takes a point $p$ on section SEC1 to section SEC2 by forward
// integration.
//
// PARAMETERS
// ==========
// mu
//    mass parameter for the RTBP
// H
//    energy value
// p
//    Initial point, 2 coordinates: p=(g,G). 
//    On return of the this function, it holds the image point on the section
//    SEC2.
// ti
//    On exit, *ti holds the integration time to reach the section SEC2.
// 
// RETURN VALUE
// ============
// Returns a non-zero error code to indicate an error and 0 to indicate
// success.
//
// NOTES
// =====
//
// CALLS TO: hinv_del, prtbp_del
//
// CALLED FROM:

int sec1sec2_del(double mu, double H, double p[2], double *ti)
{
   double x[DIM];

   x[0]=0;	// l
   x[2]=p[0];	// g
   x[3]=p[1]; 	// G

   // Compute x[1]=L by inverting the Hamiltonian.
   // Since p is on section SEC1, we need to call hinv, NOT hinv2.
   if(hinv_del(mu,H,x))
   {
      fprintf(stderr, "sec1sec2_del: error inverting the Hamiltonian\n");
      return(1);
   }

   // Compute first intersection with section SEC2 (and integration time).
   // On exit, "prtbp_del" guarantees that x is exactly on Poincare section
   // SEC2.
   if(prtbp_del(mu,SEC2,1,x,ti))
   {
      fprintf(stderr, "sec1sec2_del: error computing poincare map\n");
      return(1);
   }
   // Set the image point on section SEC2.
   p[0]=x[2]; p[1]=x[3];
   return(0);
}

int sec1sec2_del_inv(double mu, double H, double p[2], double *ti)
{
   double x[DIM];

   x[0]=0;	// l
   x[2]=p[0];	// g
   x[3]=p[1]; 	// G

   // Compute x[1]=L by inverting the Hamiltonian.
   // Since p is on section SEC1, we need to call hinv, NOT hinv2.
   if(hinv_del(mu,H,x))
   {
      fprintf(stderr, "sec1sec2_del_inv: error inverting the Hamiltonian\n");
      return(1);
   }

   // Compute first intersection with section SEC2 (and integration time).
   // On exit, "prtbp_del_inv" guarantees that x is exactly on Poincare
   // section SEC2.
   if(prtbp_del_inv(mu,SEC2,1,x,ti))
   {
      fprintf(stderr, "sec1sec2_del_inv: error computing inverse poincare map\n");
      return(1);
   }
   // Set the image point on section S2.
   p[0]=x[2]; p[1]=x[3];
   return(0);
}
