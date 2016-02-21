// ======================================================
// Poincare map of the Restricted Three Body Problem (2D)
// ======================================================
// FILE:          $RCSfile: prtbps_2d_main.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2011-02-22 11:41:08 $
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
//    - number of iterates "n"
//
// 2. For each input line:
//
//    2.1. Input data:
//    - energy value "H"
//    - initial point "x"
//
//    2.2. Compute n-th iterate of the Poincare map, $P^n(x)$.
//
//    2.3. Output final point $P^n(x)$ and integration time to stdout.

#include <stdio.h>
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE
#include <gsl/gsl_errno.h>	// gsl_set_error_handler_off
#include "prtbp_2d.h"	// prtbp2d

int main( )
{
   double mu, H, ti;
   double x[2];
   int status, n;

   // Input mass parameter, number of iterates from stdin.
   if(scanf("%le %d", &mu, &n) < 2)
   {
      perror("main: error reading input");
      exit(EXIT_FAILURE);
   }

   // Stop GSL default error handler from aborting the program
   gsl_set_error_handler_off();

   // For each energy level H in the range, do
   while(scanf("%le %le %le", &H, x, x+1)==3)
   {
      // Compute n-th iterate of 2D Poincare map, $P^n(x)$ (and integration time
      // ti).
      status=prtbp_2d(mu,H,n,x,&ti);
      if(status)
      {
	 fprintf(stderr, \
	       "main: error computing %d-th iterate of Poincare map\n",n);
	 exit(EXIT_FAILURE);
      }

      // Output final point and integration time to stdout.
      if(printf("%e %.15le %.15le %.15le\n", H, x[0], x[1], ti)<0)
      {
	 perror("main: error writting output");
	 exit(EXIT_FAILURE);
      }
   }
   exit(EXIT_SUCCESS);
}
