// ====================================================================
// Derivative of Poincare map of the Restricted Three Body Problem (2D)
// ====================================================================
// FILE:          $RCSfile: dprtbp2_2d_main.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2011-04-12 13:54:40 $
//
// PURPOSE
// =======
// Compute the derivative of 2D Poincare map $DP^n(x)$ of the RTBP.
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
//    - initial point "x"
//
// 2. Compute derivative of n-th iterate of the Poincare map, $DP^n(x)$.
// 3. Output derivative $DP^n(x)$ to stdout.

#include <stdio.h>
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE
#include <gsl/gsl_errno.h>	// gsl_set_error_handler_off
#include "dprtbp2_2d.h"	// dprtbp2d

int main( )
{
   double mu, H;
   double x[2];
   double dp[4];	// derivative of 2D Poincare map
   int status, n;

   // Input mass parameter, energy value, number of iterates, initial
   // condition, momentum estimate from stdin.
   if(scanf("%le %le %d %le %le", &mu, &H, &n, x, x+1) < 5)
   {
      perror("main: error reading input");
      exit(EXIT_FAILURE);
   }

   // Stop GSL default error handler from aborting the program
   gsl_set_error_handler_off();

   // Compute derivative of n-th iterate of 2D Poincare map, $DP^n(x)$.
   status=dprtbp2_2d(mu,H,n,x,dp);
   if(status)
   {
      fprintf(stderr, \
	    "main: error computing derivative of %d-th iterate of Poincare map\n",n);
      exit(EXIT_FAILURE);
   }

   // Output derivative to stdout.
   printf("% .15le % .15le\n", dp[0], dp[1]);
   printf("% .15le % .15le\n", dp[2], dp[3]);
   exit(EXIT_SUCCESS);
}
