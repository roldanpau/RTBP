// ===============================================
// Test Derivative of Poincare Map and its Inverse
// ===============================================
// FILE:          $RCSfile: dprtbp_test.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2012-11-26 09:48:57 $
//
// PURPOSE
// =======
// Let S be the Poincare section {y=0} in the RTBP.
// This program tests the derivative of the n-th iterate of the Poincare map,
// $DP^n(x)$, and its inverse derivative, $DP^{-n}(y)$.
// To test, we compute the derivative DP^n(x) at a point x, and the inverse
// derivative DP^{-n}(y) at y=P^n(x), and then check that the product
//    DP^{-n}(y)*DP^n(x) = Id
// is equal to the identity.
//
// NOTES
// =====
//
// OVERALL METHOD:
//
// 1. Input mass parameter, number of iterates "n", and argument to the
//    derivative "x" from stdin.
// 2. Compute both direct and inverse derivative, and multiply.
// 3. Output the product DP^{-n}(y)*DP^n(x) to stdout, which should be the
// identity.

#include <stdio.h>
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE
#include <gsl/gsl_blas.h>	// gsl_blas_gemm
#include <gsl/gsl_errno.h>	// gsl_set_error_handler_off
#include <rtbp.h>	// DIM
#include <prtbp.h>	// prtbp
#include "dprtbp.h"	// dprtbp

int main( )
{
   double mu;
   double x[DIM], y[DIM];
   double dp[DIMV], dp_inv[DIMV];	// derivative of poincare map
   double id[DIMV];			// product matrix id=dp*dp_inv
   double ti;				// integration time
   int status, n;

   // Input mass parameter, number of iterates and initial condition from stdin.
   if(scanf("%le %d %le %le %le %le", &mu, &n, x, x+1, x+2, x+3) < 6)
   {
      perror("main: error reading input");
      exit(EXIT_FAILURE);
   }

   // Stop GSL default error handler from aborting the program
   gsl_set_error_handler_off();

   // Compute derivative of n-th iterate of Poincare map, $DP^n(x)$.
   status=dprtbp(mu,SEC1,n,x,dp);
   if(status)
   {
      fprintf(stderr, \
	    "main: error computing derivative of %d-th iterate of Poincare map\n",n);
      exit(EXIT_FAILURE);
   }

   // Compute derivative of n-th iterate of inverse Poincare map, $DP^{-n}(y)$.

   // y=P^{n}(x)
   y[0]=x[0]; y[1]=x[1]; y[2]=x[2]; y[3]=x[3];
   if(prtbp(mu,SEC1,n,y,&ti))
   {
      fprintf(stderr, "main: error computing poincare map\n");
      exit(EXIT_FAILURE);
   }

   status=dprtbp_inv(mu,SEC1,n,y,dp_inv);
   if(status)
   {
      fprintf(stderr, \
	    "main: error computing derivative of %d-th iterate of inverse Poincare map\n",n);
      exit(EXIT_FAILURE);
   }

// 3. Output the product DP^{-n}(y)*DP^n(x) to stdout, which should be the
// identity.
   gsl_matrix_view A=gsl_matrix_view_array(dp,DIM,DIM);
   gsl_matrix_view B=gsl_matrix_view_array(dp_inv,DIM,DIM);
   gsl_matrix_view C=gsl_matrix_view_array(id,DIM,DIM);

   // Compute A= B C
   gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                       1.0, &A.matrix, &B.matrix,
                       0.0, &C.matrix);

   // Output product to stdout.
   printf("% .15le % .15le % .15le % .15le\n", id[ 0], id[ 1], id[ 2], id[ 3]);
   printf("% .15le % .15le % .15le % .15le\n", id[ 4], id[ 5], id[ 6], id[ 7]);
   printf("% .15le % .15le % .15le % .15le\n", id[ 8], id[ 9], id[10], id[11]);
   printf("% .15le % .15le % .15le % .15le\n", id[12], id[13], id[14], id[15]);
   exit(EXIT_SUCCESS);
}
