// ==============================================================
// Deivative of Poincare map of the Restricted Three Body Problem
// ==============================================================
// FILE:          $RCSfile: dprtbp_main.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2012-11-26 09:48:57 $
//
// PURPOSE
// =======
// Let S be the Poincare section {y=0} in the RTBP.
// This program computes the derivative of the n-th iterate of the Poincare
// map, $DP^n(x)$.
//
// NOTES
// =====
// We do not check (and thus we do not impose) that the initial point "x" is
// on the Poincare section.
//
// OVERALL METHOD:
//
// 1. Input mass parameter, Poincare section "sec", number of iterates "n",
//    and argument to the derivative "x" from stdin.
// 2. Compute derivative of n-th iterate of the Poincare map, $DP^n(x)$.
// 3. Output derivative $DP^n(x)$ to stdout.

#include <stdio.h>
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE
#include <gsl/gsl_errno.h>	// gsl_set_error_handler_off
#include <rtbp.h>	// DIM
#include <prtbp.h>	// section_t
#include "dprtbp.h"	// dprtbp

int main( )
{
   double mu;
   section_t sec;	// type of Poincare section
   double x[DIM];
   double dp[DIMV];	// derivative of poincare map
   int status, n;

   // auxiliary variables
   char section_str[10];        // holds input string "SEC1", "SEC2" etc

   // Input mass parameter, type of section, number of iterates and initial
   // condition from stdin.
   if(scanf("%le %s %d %le %le %le %le", &mu, section_str, &n, x, x+1, x+2, x+3) < 7)
   {
      perror("main: error reading input");
      exit(EXIT_FAILURE);
   }

   // Stop GSL default error handler from aborting the program
   gsl_set_error_handler_off();

   // Compute derivative of n-th iterate of Poincare map, $DP^n(x)$.
   status=dprtbp(mu,sec,n,x,dp);
   if(status)
   {
      fprintf(stderr, \
	    "main: error computing derivative of %d-th iterate of Poincare map\n",n);
      exit(EXIT_FAILURE);
   }

   // Output derivative to stdout.
   printf("% .15le % .15le % .15le % .15le\n", dp[ 0], dp[ 1], dp[ 2], dp[ 3]);
   printf("% .15le % .15le % .15le % .15le\n", dp[ 4], dp[ 5], dp[ 6], dp[ 7]);
   printf("% .15le % .15le % .15le % .15le\n", dp[ 8], dp[ 9], dp[10], dp[11]);
   printf("% .15le % .15le % .15le % .15le\n", dp[12], dp[13], dp[14], dp[15]);
   exit(EXIT_SUCCESS);
}
