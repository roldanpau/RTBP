// ======================================================
// Poincare map of the Restricted Three Body Problem (2D)
// ======================================================
// FILE:          $RCSfile: prtbp2_2d_main.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2011-03-09 15:36:59 $
//
// PURPOSE
// =======
// Compute the 2D Poincare map $P(x)$ of the RTBP.
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
//
// 2. For each input line:
//
//    2.1. Input data:
//    - initial point "x"
//
//    2.2. Compute n-th iterate of the Poincare map, $P^n(x)$.
//
//    2.3. Output final point $P^n(x)$ and integration time to stdout.

#include <stdio.h>
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE
#include <gsl/gsl_errno.h>	// gsl_set_error_handler_off
#include "prtbp2_2d.h"	// prtbp2d

int main( )
{
   double mu, H, ti;
   double x[2];
   int status, n;

   // Input mass parameter, energy value, number of iterates from stdin.
   if(scanf("%le %le %d", &mu, &H, &n) < 3)
   {
      perror("main: error reading input");
      exit(EXIT_FAILURE);
   }

   // Stop GSL default error handler from aborting the program
   gsl_set_error_handler_off();

   // For each initial condition, do
   while(scanf("%le %le", x, x+1)==2)
   {
   // Compute n-th iterate of 2D Poincare map, $P^n(x)$ (and integration time
   // ti).
   status=prtbp2_2d(mu,H,n,x,&ti);
   if(status)
   {
      fprintf(stderr, \
	    "main: error computing %d-th iterate of Poincare map\n",n);
      exit(EXIT_FAILURE);
   }

   // Output final point and integration time to stdout.
   if(printf("%.15le %.15le %.15le\n", x[0], x[1], ti)<0)
   {
      perror("main: error writting output");
      exit(EXIT_FAILURE);
   }
   }
   exit(EXIT_SUCCESS);
}
