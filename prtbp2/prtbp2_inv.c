// =========================================================
// Inverse Poincare map of the Restricted Three Body Problem
// =========================================================
// FILE:          $RCSfile: prtbp2_inv_main.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2011-03-09 15:36:59 $
//
// PURPOSE
// =======
// Let S be the Poincare section {y=0} in the RTBP.
// This program computes the inverse Poincare map, $P^{-1}(x)$.
//
// NOTES
// =====
// We do not check (and thus we do not impose) that the initial point "x" is
// on the Poincare section.
//
// OVERALL METHOD:
//
// 1. Input mass parameter, number of iterates "n" and initial point "x" from stdin.
// 2. Compute inverse Poincare map, $P^{-1}(x)$.
// 3. Output final point $P^{-1}(x)$ and integration time to stdout.

#include <stdio.h>
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE
#include <gsl/gsl_errno.h>	// gsl_set_error_handler_off
#include <rtbp.h>	// DIM
#include "prtbp2.h"	// prtbp2_inv

int main( )
{
   double mu;
   double ti;	// integration time
   double x[DIM];
   int status, n;

   // Input mass parameter, number of iterates and initial condition from stdin.
   if(scanf("%le %d %le %le %le %le", &mu, &n, x, x+1, x+2, x+3) < 6)
   {
      perror("main: error reading input");
      exit(EXIT_FAILURE);
   }

   // Stop GSL default error handler from aborting the program
   gsl_set_error_handler_off();

   // Compute inverse Poincare map, $P^{-1}(x)$.
   status=prtbp2_inv(mu,n,x,&ti);
   if(status)
   {
      fprintf(stderr, "main: error computing inverse Poincare map\n");
      exit(EXIT_FAILURE);
   }

   // Output final point and integration time to stdout.
   status = printf("%.15le %.15le %.15le %.15le %.15le\n", \
	 x[0], x[1], x[2], x[3], ti);
   if(status<0)
   {
      perror("main: error writting output");
      exit(EXIT_FAILURE);
   }
   exit(EXIT_SUCCESS);
}
