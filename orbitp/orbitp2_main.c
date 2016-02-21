// =============================
// Orbit of Poincare map in RTBP
// =============================
// FILE:          $RCSfile: orbitp2_main.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2012-05-02 09:43:04 $
//
// PURPOSE
// =======
// Consider the Restricted Three Body Problem. Assume that energy is fixed:
// $H(x,y,p_x,p_y)=\bar H$.
// Let $\Sigma^- = {y=0, p_y<0}$ be a Poincare section, and $P(x,p_x)$ be the
// associated 2D Poincare map.
// Let $\tilde P = P^k$ be the iterated Poincare map.
// Given a point $p=(x,p_x)$, this function computes the orbit 
// ${\tilde P^i(p)}$ for $i=1,2,\dots,n$.
//
// OVERALL METHOD
// ==============
//
// 1. Input parameters from stdin:
// 
//    - mass parameter 
//    - energy value "H"
//    - number of iterates "k" of iterated map
//    - initial point "x"
//    - number of desired points in the orbit
//
// 2. Compute orbit of $x$ under the Poincare map: ${\tilde P^i(x)}$.
// 3. Output orbit to stdout.

#include <stdio.h>
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE
#include <gsl/gsl_errno.h>	// gsl_set_error_handler_off
#include "orbitp2.h"	// orbitp2, print_orbit

int main( )
{
   double mu, H;
   int k;
   double x[2];
   int npt;
   double *orbit;

   // Input mass parameter, energy value, number of iterates, initial
   // condition, and number of points in the orbit from stdin.
   if(scanf("%le %le %d %le %le %d", &mu, &H, &k, x, x+1, &npt) < 6)
   {
      perror("main: error reading input");
      exit(EXIT_FAILURE);
   }

   orbit = (double *)malloc(npt*2*sizeof(double));

   // Stop GSL default error handler from aborting the program
   gsl_set_error_handler_off();

   orbitp2(mu,H,k,x,npt,orbit);

   print_orbit(npt,orbit);
   exit(EXIT_SUCCESS);
}
