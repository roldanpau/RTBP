// ===========================================================
// Hyperbolic Splitting Associated to a Family of Fixed Points
// ===========================================================
// FILE:          $RCSfile: hypers_main.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2013-02-20 09:59:41 $
//
// PURPOSE
// =======
// Consider a 2D map $F$ and its derivative $DF$.
// Let $\Lambda$ be a family of hyperbolic fixed points for $F$.
// For every fixed point $x \in \Lambda$, compute the hyperbolic splitting,
// i.e. its unstable/stable eigenvalues and eigenvectors.
//
// NOTES
// =====
//
// OVERALL METHOD:
//
// 1. Input parameters from stdin:
// 
//    - mass parameter 
//    - Poincare section "sec"
//    - number of iterates "n"
//
// 2. Each line in stdin corresponds to a fixed point $x$. 
//    For each energy level H in the range, do
//    2.1. Input one line in stdin: energy value "H", fixed point "x".
//    2.2. Compute hyperbolic splitting associated to the fixed point $x$.
//    2.3. Output eigenvalues and eigenvectors of $x$ to stdout.

#include <stdio.h>
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE
#include <string.h>	// strcmp
#include <gsl/gsl_errno.h>	// gsl_set_error_handler_off
#include <prtbp.h>	// section_t
#include "hyper.h"	// hyper

int main( )
{
   double mu, H;
   section_t sec;
   double x[2];
   double eval[2], evec[4];	// eigenvalues/eigenvectors
   int status, n;

   // auxiliary variables
   char section_str[10];        // holds input string "SEC1", "SEC2" etc

   // Input mass parameter, Poincare section, number of iterates from stdin.
   if(scanf("%le %s %d", &mu, section_str, &n) < 3)
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

   // For each energy level H in the range, do
   while(scanf("%le %le %le", &H, x, x+1)==3)
   {
      // Compute derivative of n-th iterate of 2D Poincare map, $DP^n(x)$.
      status=hyper(mu,sec,H,n,x,eval,evec);
      if(status)
      {
	 fprintf(stderr, \
	       "main: error computing hyperbolic splitting\n");
	 exit(EXIT_FAILURE);
      }

      // Output eigenvalues, eigenvectors to stdout.
      printf("%.15le ", H);
      printf("%.15le %.15le ", eval[0], eval[1]);	// lambda_u, _s
      printf("%.15le %.15le ", evec[0], evec[1]);	// v_u[2]
      printf("%.15le %.15le\n", evec[2], evec[3]);	// v_s[2]
   }
   exit(EXIT_SUCCESS);
}
