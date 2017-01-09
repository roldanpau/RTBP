// =========================================================
// Inverse Poincare map of the Restricted Three Body Problem
// =========================================================
// FILE:          $RCSfile: prtbp_inv_main.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2012-11-26 10:39:05 $
//
// PURPOSE
// =======
// Let "sec" be a Poincare section in the RTBP, where
//      - sec = SEC1 means section {y=0, p_y>0}
//      - sec = SEC2 means section {y=0, p_y<0}.
// This program computes the n-th iterate of the inverse Poincare map,
// $P^{-n}(x)$.
//
// NOTES
// =====
// We do not check (and thus we do not impose) that the initial point "x" is
// on the Poincare section.
//
// OVERALL METHOD:
//
// 1. Input mass parameter, type of section "sec", number of iterates "n" and
// initial point "x" from stdin.
// 2. Compute n-th iterate of the inverse Poincare map, $P^{-n}(x)$.
// 3. Output final point $P^{-n}(x)$ and integration time to stdout.

#include <stdio.h>
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE
#include <string.h>     // strcmp

#include <gsl/gsl_errno.h>	// gsl_set_error_handler_off

#include <rtbp.h>	// DIM
#include "prtbp.h"	// section_t, prtbp_inv

int main( )
{
   double mu;
   section_t sec;
   double ti;		// integration time
   double x[DIM];
   int status, n;

   // auxiliary variables
   char section_str[10];        // holds input string "SEC1", "SEC2" etc

   // Input mass parameter, section, number of iterates and initial condition from stdin.
   if(scanf("%le %s %d %le %le %le %le", &mu, section_str, &n, x, x+1, x+2, x+3) < 7)
   {
      perror("main: error reading input");
      exit(EXIT_FAILURE);
   }

   if (strcmp(section_str,"SEC1") == 0)
      sec = SEC1;
   else if (strcmp(section_str,"SEC2") == 0)
      sec = SEC2;
   else
   {
      perror("main: error reading section string");
      exit(EXIT_FAILURE);
   }

   // Stop GSL default error handler from aborting the program
   gsl_set_error_handler_off();

   // Compute n-th iterate of inverse Poincare map, $P^{-n}(x)$.
   status=prtbp_inv(mu,sec, n,x,&ti);
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
