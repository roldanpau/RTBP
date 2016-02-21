// ====================================================================
// Derivative of Poincare map of the Restricted Three Body Problem (2D)
// ====================================================================
// FILE:          $RCSfile: dprtbpdel_2d_main.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2012-05-31 06:52:10 $
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
//    - type of poincare section "sec"
//    - energy value "H"
//    - number of iterates "n"
//    - initial point "x"
//
// 2. Compute derivative of n-th iterate of the Poincare map, $DP^n(x)$.
// 3. Output derivative $DP^n(x)$ to stdout.

#include <stdio.h>
#include <stdlib.h>		// EXIT_SUCCESS, EXIT_FAILURE
#include <gsl/gsl_errno.h>	// gsl_set_error_handler_off
#include <prtbpdel.h>		// section_t
#include "dprtbpdel_2d.h"	// dprtbp_del_2d

int main( )
{
   double mu, H;
   section_t sec;	// type of Poincare section
   double x[2];
   double dp[4];	// derivative of 2D Poincare map
   int status, n;

   // auxiliary variables
   char section_str[10];        // holds input string "SEC1", "SEC2" etc

   // Input mass parameter, type of section, energy value, number of
   // iterates, initial condition, momentum estimate from stdin.
   if(scanf("%le %s %le %d %le %le", &mu, section_str, &H, &n, x, x+1) < 6)
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

   // Compute derivative of n-th iterate of 2D Poincare map, $DP^n(x)$.
   status=dprtbp_del_2d(mu,sec,H,n,x,dp);
   if(status)
   {
      fprintf(stderr, \
	    "main: error computing derivative of %d-th iterate of Poincare map\n",n);
      exit(EXIT_FAILURE);
   }

   // Output derivative to stdout.
   printf("% .10le % .10le\n", dp[0], dp[1]);
   printf("% .10le % .10le\n", dp[2], dp[3]);
   exit(EXIT_SUCCESS);
}
