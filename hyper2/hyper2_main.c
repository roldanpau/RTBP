// ====================================================================
// Hyperbolic Splitting (Eigenvalues/vects) Associated to a Fixed Point
// ====================================================================
// FILE:          $RCSfile: hyper2_main.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2011-03-17 16:12:11 $
//
// PURPOSE
// =======
// Consider a 2D map $F$ and its derivative $DF$.
// Given a hyperbolic fixed point $x$ of $F$, compute its hyperbolic
// splitting, i.e. its unstable/stable eigenvalues and eigenvectors. 
//
// NOTES
// =====
//
// OVERALL METHOD:
//
// 1. Input parameters from stdin:
// 
//    - mass parameter 
//    - energy value "H"
//    - number of iterates "n"
//    - fixed point "x"
//
// 2. Compute hyperbolic splitting associated to the fixed point $$.
// 3. Output eigenvalues and eigenvectors to stdout.

#include <stdio.h>
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE
#include <gsl/gsl_errno.h>	// gsl_set_error_handler_off
#include "hyper2.h"	// hyper2

int main( )
{
   double mu, H;
   double x[2];
   double eval[2], evec[4];	// eigenvalues/eigenvectors
   int status, n;

   // Input mass parameter, energy value, number of iterates, fixed point
   // from stdin.
   if(scanf("%le %le %d %le %le", &mu, &H, &n, x, x+1) < 5)
   {
      perror("main: error reading input");
      exit(EXIT_FAILURE);
   }

   // Stop GSL default error handler from aborting the program
   gsl_set_error_handler_off();

   // Compute derivative of n-th iterate of 2D Poincare map, $DP^n(x)$.
   status=hyper2(mu,H,n,x,eval,evec);
   if(status)
   {
      fprintf(stderr, \
	    "main: error computing hyperbolic splitting\n",n);
      exit(EXIT_FAILURE);
   }

   // Output eigenvalues to stdout.
   printf("rho_u = % .15le\n", eval[0]);
   printf("rho_s = % .15le\n", eval[1]);
   printf("v_u = (% .15le, % .15le)\n", evec[0], evec[1]);
   printf("v_s = (% .15le, % .15le)\n", evec[2], evec[3]);
   exit(EXIT_SUCCESS);
}
