// ===========================================================
// Hyperbolic Splitting Associated to a Family of Fixed Points
// ===========================================================
// FILE:          $RCSfile: hyper2s_main.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2011-05-05 13:22:19 $
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
// There are two possible orientations for the eigenvectors. 
// We choose the following convention (because we are thinking of the
// manifolds corresponding to the "inner" separatrix of the pendulum
// associated to the fixed points $p_2,p_3$ on the Poincare section $\Sigma'$):
// - For the unstable eigenvector, we choose p_x>0.
// - For the stable eigenvector, we choose p_x<0.
//
// OVERALL METHOD:
//
// 1. Input parameters from stdin:
// 
//    - mass parameter 
//    - number of iterates "n"
//
// 2. Each line in stdin corresponds to a fixed point $x$. 
//    For each energy level H in the range, do
//    2.1. Input one line in stdin: energy value "H", fixed point "x".
//    2.2. Compute hyperbolic splitting associated to the fixed point $x$.
//    2.3. Output eigenvalues and eigenvectors of $x$ to stdout.

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
      // Compute derivative of n-th iterate of 2D Poincare map, $DP^n(x)$.
      status=hyper2(mu,H,n,x,eval,evec);
      if(status)
      {
	 fprintf(stderr, \
	       "main: error computing hyperbolic splitting\n",n);
	 exit(EXIT_FAILURE);
      }

      // - For the unstable eigenvector, we choose p_x>0.
      if (evec[1]<0)
      {
	 evec[0] = -evec[0];
	 evec[1] = -evec[1];
      }
      // - For the stable eigenvector, we choose p_x<0.
      if (evec[3]>0)
      {
	 evec[2] = -evec[2];
	 evec[3] = -evec[3];
      }

      // Output eigenvalues, eigenvectors to stdout.
      printf("%.15le ", H);
      printf("%.15le %.15le ", eval[0], eval[1]);	// lambda_u, _s
      printf("%.15le %.15le ", evec[0], evec[1]);	// v_u[2]
      printf("%.15le %.15le\n", evec[2], evec[3]);	// v_s[2]
   }
   exit(EXIT_SUCCESS);
}
