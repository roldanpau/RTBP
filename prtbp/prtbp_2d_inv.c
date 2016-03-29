// ==============================================================
// Inverse Poincare map of the Restricted Three Body Problem (2D)
// ==============================================================
// FILE:          $RCSfile: prtbp_2d_inv_main.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2012-11-26 10:39:05 $
//
// PURPOSE
// =======
// Compute the inverse 2D Poincare map $P^{-1}(x)$ of the RTBP.
//
// NOTES
// =====
//
// OVERALL METHOD:
//
// 1. Input parameters from stdin:
// 
//    - mass parameter 
//    - type of Poincare section "sec"
//    - energy value "H"
//    - number of iterates "n"
//
// 2. For each input line:
//
//    2.1. Input data:
//    - initial point "x"
//
//    2.2. Compute n-th iterate of the inverse Poincare map, $P^{-n}(x)$.
//
//    2.3. Output final point $P^{-n}(x)$ and integration time to stdout.

#include <stdio.h>
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE
#include <gsl/gsl_errno.h>	// gsl_set_error_handler_off
#include "prtbp.h"	// section_t
#include "prtbp_2d.h"	// prtbp_2d

int main( )
{
   double mu, H, ti;
   section_t sec;
   double x[2];
   int status, n;

   // auxiliary variables
   char section_str[10];        // holds input string "SEC1", "SEC2" etc

   // Input mass parameter, energy value, number of iterates from stdin.
   if(scanf("%le %s %le %d", &mu, section_str, &H, &n) < 4)
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

   // For each initial condition, do
   while(scanf("%le %le", x, x+1)==2)
   {
   // Compute n-th iterate of 2D Poincare map, $P^n(x)$ (and integration time
   // ti).
   status=prtbp_2d_inv(mu,sec,H,n,x,&ti);
   if(status)
   {
      fprintf(stderr, \
	    "main: error computing %d-th iterate of inverse Poincare map\n",n);
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
