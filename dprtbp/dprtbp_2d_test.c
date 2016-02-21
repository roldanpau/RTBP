// =========================================================================
// Test Derivative of Poincare map of the Restricted Three Body Problem (2D)
// =========================================================================
// FILE:          $RCSfile: dprtbp_2d_test.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2012-11-26 09:48:57 $
//
// PURPOSE
// =======
// Test the derivative of 2D Poincare map $DP^n(x)$ of the RTBP by comparing
// it to the derivative obtained using numerical differentiation.
//
// Test also the derivative of inverse 2D Poicare map by verifying the
// product 
//    DP^n(x) DP^{-n}(y) = Id,
// where y=P^n(x).
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
//    - initial point "x"
//
// 2. Compute derivative of n-th iterate of the Poincare map, $DP^n(x)$.
// 3. Compute the same derivative using numerical differentiation.
// 4. Output the difference of both matrices to stdout.
// 5. Compute $DP^{-n}(y)$, where y=P^n(x).
// 6. Output the product DP^n(x) DP^{-n}(y) to stdout (it should be the
// identity).

#include <stdio.h>
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE
#include <gsl/gsl_errno.h>	// gsl_set_error_handler_off
#include <gsl/gsl_deriv.h>	// gsl_deriv_central
#include <gsl/gsl_blas.h>	// gsl_blas_gemm
#include <prtbp_2d.h>	// prtbp_2d, prtbp_2d_inv
#include "dprtbp_2d.h"	// dprtbp_2d

struct dparams
{
   double mu;
   double H;
   double arg;	// fixed value of one argument of Poincare map
   int cuts;	// number of iterates of poincare map
};

// First component of Poincare map as a function of the first variable (x)

double p11(double x, void *params)
{
   double mu=((struct dparams *)params)->mu;
   double H=((struct dparams *)params)->H;
   int n=((struct dparams *)params)->cuts;
   double pt[2];

   pt[0]=x;
   pt[1]=((struct dparams *)params)->arg;	// p_x

   // this is not used, but it's needed in call to prtbp_2d
   double ti;   // integration time

   // Compute the 2D poincare map $P^(pt)$
   if(prtbp_2d(mu,SEC1,H,n,pt,&ti))
   {
      fprintf(stderr, "dist: error computing 2D poincare map\n");
      exit(1);
   }
   return(pt[0]);
}

// First component of Poincare map as a function of the second variable (px)
double p12(double px, void *params)
{
   double mu=((struct dparams *)params)->mu;
   double H=((struct dparams *)params)->H;
   int n=((struct dparams *)params)->cuts;
   double pt[2];

   pt[0]=((struct dparams *)params)->arg;	// x
   pt[1]=px;

   // this is not used, but it's needed in call to prtbp_2d
   double ti;   // integration time

   // Compute the 2D poincare map $P^(pt)$
   if(prtbp_2d(mu,SEC1,H,n,pt,&ti))
   {
      fprintf(stderr, "dist: error computing 2D poincare map\n");
      exit(1);
   }
   return(pt[0]);
}

// Second component of Poincare map as a function of the first variable (x)
double p21(double x, void *params)
{
   double mu=((struct dparams *)params)->mu;
   double H=((struct dparams *)params)->H;
   int n=((struct dparams *)params)->cuts;
   double pt[2];

   pt[0]=x;
   pt[1]=((struct dparams *)params)->arg;	// p_x

   // this is not used, but it's needed in call to prtbp_2d
   double ti;   // integration time

   // Compute the 2D poincare map $P^(pt)$
   if(prtbp_2d(mu,SEC1,H,n,pt,&ti))
   {
      fprintf(stderr, "dist: error computing 2D poincare map\n");
      exit(1);
   }
   return(pt[1]);
}

// Second component of Poincare map as a function of the second variable (px)
double p22(double px, void *params)
{
   double mu=((struct dparams *)params)->mu;
   double H=((struct dparams *)params)->H;
   int n=((struct dparams *)params)->cuts;
   double pt[2];

   pt[0]=((struct dparams *)params)->arg;	// x
   pt[1]=px;

   // this is not used, but it's needed in call to prtbp_2d
   double ti;   // integration time

   // Compute the 2D poincare map $P^(pt)$
   if(prtbp_2d(mu,SEC1,H,n,pt,&ti))
   {
      fprintf(stderr, "dist: error computing 2D poincare map\n");
      exit(1);
   }
   return(pt[1]);
}

// First component of Poincare map as a function of the first variable (x)

double p11_inv(double x, void *params)
{
   double mu=((struct dparams *)params)->mu;
   double H=((struct dparams *)params)->H;
   int n=((struct dparams *)params)->cuts;
   double pt[2];

   pt[0]=x;
   pt[1]=((struct dparams *)params)->arg;	// p_x

   // this is not used, but it's needed in call to prtbp_2d
   double ti;   // integration time

   // Compute the 2D poincare map $P^(pt)$
   if(prtbp_2d_inv(mu,SEC1,H,n,pt,&ti))
   {
      fprintf(stderr, "dist: error computing 2D poincare map\n");
      exit(1);
   }
   return(pt[0]);
}

// First component of Poincare map as a function of the second variable (px)
double p12_inv(double px, void *params)
{
   double mu=((struct dparams *)params)->mu;
   double H=((struct dparams *)params)->H;
   int n=((struct dparams *)params)->cuts;
   double pt[2];

   pt[0]=((struct dparams *)params)->arg;	// x
   pt[1]=px;

   // this is not used, but it's needed in call to prtbp_2d
   double ti;   // integration time

   // Compute the 2D poincare map $P^(pt)$
   if(prtbp_2d_inv(mu,SEC1,H,n,pt,&ti))
   {
      fprintf(stderr, "dist: error computing 2D poincare map\n");
      exit(1);
   }
   return(pt[0]);
}

// Second component of Poincare map as a function of the first variable (x)
double p21_inv(double x, void *params)
{
   double mu=((struct dparams *)params)->mu;
   double H=((struct dparams *)params)->H;
   int n=((struct dparams *)params)->cuts;
   double pt[2];

   pt[0]=x;
   pt[1]=((struct dparams *)params)->arg;	// p_x

   // this is not used, but it's needed in call to prtbp_2d
   double ti;   // integration time

   // Compute the 2D poincare map $P^(pt)$
   if(prtbp_2d_inv(mu,SEC1,H,n,pt,&ti))
   {
      fprintf(stderr, "dist: error computing 2D poincare map\n");
      exit(1);
   }
   return(pt[1]);
}

// Second component of Poincare map as a function of the second variable (px)
double p22_inv(double px, void *params)
{
   double mu=((struct dparams *)params)->mu;
   double H=((struct dparams *)params)->H;
   int n=((struct dparams *)params)->cuts;
   double pt[2];

   pt[0]=((struct dparams *)params)->arg;	// x
   pt[1]=px;

   // this is not used, but it's needed in call to prtbp_2d
   double ti;   // integration time

   // Compute the 2D poincare map $P^(pt)$
   if(prtbp_2d_inv(mu,SEC1,H,n,pt,&ti))
   {
      fprintf(stderr, "dist: error computing 2D poincare map\n");
      exit(1);
   }
   return(pt[1]);
}

int main( )
{
   double mu, H, ti;
   double x[2], y[2];
   double dp[4], dp_inv[4];	// derivative of 2D Poincare map
   double dpn[4], dpn_inv[4];	// derivative of 2D Poincare map (numerical)
   double id[4];		// product matrix id=dp*dp_inv
   double abserr[4];		// error of numerical derivative
   int status, n;
   struct dparams params;

   gsl_function P11, P12, P21, P22;

   // Input mass parameter, energy value, number of iterates, initial
   // condition from stdin.
   if(scanf("%le %le %d %le %le", &mu, &H, &n, x, x+1) < 5)
   {
      perror("main: error reading input");
      exit(EXIT_FAILURE);
   }

   // Stop GSL default error handler from aborting the program
   gsl_set_error_handler_off();

   // Compute derivative of n-th iterate of 2D Poincare map, $DP^n(x)$.
   status=dprtbp_2d(mu,SEC1,H,n,x,dp);
   if(status)
   {
      fprintf(stderr, \
	    "main: error computing derivative of %d-th iterate of Poincare map\n",n);
      exit(EXIT_FAILURE);
   }

   // Compute the same derivative using numerical differentiation.
   params.mu = mu;
   params.H = H;
   params.cuts = n;

   params.arg = x[1];
   P11.function = &p11;
   P11.params = &params;
   gsl_deriv_central(&P11, x[0], 1.e-08, dpn, abserr);

   params.arg = x[0];
   P12.function = &p12;
   P12.params = &params;
   gsl_deriv_central(&P12, x[1], 1.e-08, dpn+1, abserr+1);

   params.arg = x[1];
   P21.function = &p21;
   P21.params = &params;
   gsl_deriv_central(&P21, x[0], 1.e-08, dpn+2, abserr+2);

   params.arg = x[0];
   P22.function = &p22;
   P22.params = &params;
   gsl_deriv_central(&P22, x[1], 1.e-08, dpn+3, abserr+3);

   printf("Estimate of absolute error in numerical derivative:\n");
   printf("% .10le % .10le\n", abserr[0], abserr[1]);
   printf("% .10le % .10le\n", abserr[2], abserr[3]);

   // Output the difference of both matrices to stdout.
   printf("Discrepancy between analytical and numerical derivative:\n");
   printf("% .10le % .10le\n", dp[0]-dpn[0], dp[1]-dpn[1]);
   printf("% .10le % .10le\n", dp[2]-dpn[2], dp[3]-dpn[3]);
   printf("\n");

   // y=P^{n}(x)
   y[0]=x[0]; y[1]=x[1];
   if(prtbp_2d(mu,SEC1,H,n,y,&ti))
   {
      fprintf(stderr, "dist: error computing 2D poincare map\n");
      exit(EXIT_FAILURE);
   }

   // Compute derivative of n-th iterate of inverse Poincare map, $DP^{-n}(y)$.
   status=dprtbp_2d_inv(mu,SEC1,H,n,y,dp_inv);
   if(status)
   {
      fprintf(stderr, \
	    "main: error computing derivative of %d-th iterate of inverse Poincare map\n",n);
      exit(EXIT_FAILURE);
   }

   // Compute the same derivative using numerical differentiation.
   params.mu = mu;
   params.H = H;
   params.cuts = n;

   params.arg = y[1];
   P11.function = &p11_inv;
   P11.params = &params;
   gsl_deriv_central(&P11, y[0], 1.e-08, dpn_inv, abserr);

   params.arg = y[0];
   P12.function = &p12_inv;
   P12.params = &params;
   gsl_deriv_central(&P12, y[1], 1.e-08, dpn_inv+1, abserr+1);

   params.arg = y[1];
   P21.function = &p21_inv;
   P21.params = &params;
   gsl_deriv_central(&P21, y[0], 1.e-08, dpn_inv+2, abserr+2);

   params.arg = y[0];
   P22.function = &p22_inv;
   P22.params = &params;
   gsl_deriv_central(&P22, y[1], 1.e-08, dpn_inv+3, abserr+3);

   printf("Estimate of absolute error in numerical derivative:\n");
   printf("% .10le % .10le\n", abserr[0], abserr[1]);
   printf("% .10le % .10le\n", abserr[2], abserr[3]);

   // Output the difference of both matrices to stdout.
   printf("Discrepancy between analytical and numerical derivative:\n");
   printf("% .10le % .10le\n", dp_inv[0]-dpn_inv[0], dp_inv[1]-dpn_inv[1]);
   printf("% .10le % .10le\n", dp_inv[2]-dpn_inv[2], dp_inv[3]-dpn_inv[3]);
   printf("\n");

// Output the product DP^{-n}(y)*DP^n(x) to stdout, which should be the
// identity.
   gsl_matrix_view A=gsl_matrix_view_array(dp,2,2);
   gsl_matrix_view B=gsl_matrix_view_array(dp_inv,2,2);
   gsl_matrix_view C=gsl_matrix_view_array(id,2,2);

   // Compute A= B C
   gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                       1.0, &A.matrix, &B.matrix,
                       0.0, &C.matrix);

   // Output product to stdout.
   printf("% .15le % .15le\n", id[ 0], id[ 1]);
   printf("% .15le % .15le\n", id[ 2], id[ 3]);

   exit(EXIT_SUCCESS);
}
