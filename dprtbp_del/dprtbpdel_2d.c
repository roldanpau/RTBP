// ====================================================================
// Derivative of Poincare map of the Restricted Three Body Problem (2D)
// ====================================================================
// FILE:          $RCSfile: dprtbpdel_2d.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2012-05-31 06:52:10 $
//
// FUNCTIONS
// =========
//
// dprtbp_2d
// ---------
// Compute the derivative $DP$ of the 2D Poincare map.
//
// dprtbp_2d_inv
// -------------
// Compute the derivative $DP^{-1}$ of the inverse 2D Poincare map.
//
// set_dprtbp_2d
// -------------
// Compute the 2D derivative $DP(x,p_x)$ from the 4D variationals
// $DP(x,y,p_x,p_y)$, and from the vectorfield of the RTBP evaluated at the
// point $p=(x,0,p_x,p_y)$ and $\phi_T(p)$.

#include <stdio.h>		// fprintf
#include <stdlib.h>     	// exit, EXIT_FAILURE
#include <prtbpdel.h>		// section_t
#include <prtbpdel_2d.h>	// prtbp_del_2d
#include <gsl/gsl_deriv.h>	// gsl_deriv_central

const double STEP_SIZE = 1.e-8;	// for central differences algorithm

struct dparams
{
   double mu;
   section_t sec;	// type of Poincare section 
   double H;
   double arg;	// fixed value of one argument of Poincare map
   int cuts;	// number of iterates of poincare map
};

// First component of Poincare map as a function of the first variable (g)

double p11(double g, void *params)
{
   double mu=((struct dparams *)params)->mu;
   section_t sec=((struct dparams *)params)->sec;
   double H=((struct dparams *)params)->H;
   int n=((struct dparams *)params)->cuts;
   double pt[2];

   pt[0]=g;
   pt[1]=((struct dparams *)params)->arg;	// G

   // this is not used, but it's needed in call to prtbp_del_2d
   double ti;   // integration time

   // Compute the 2D poincare map $P^(pt)$
   if(prtbp_del_2d(mu,sec,H,n,pt,&ti))
   {
      fprintf(stderr, "p11: error computing 2D poincare map\n");
      exit(EXIT_FAILURE);
   }
   return(pt[0]);
}

// First component of Poincare map as a function of the second variable (G)
double p12(double G, void *params)
{
   double mu=((struct dparams *)params)->mu;
   section_t sec=((struct dparams *)params)->sec;
   double H=((struct dparams *)params)->H;
   int n=((struct dparams *)params)->cuts;
   double pt[2];

   pt[0]=((struct dparams *)params)->arg;	// g
   pt[1]=G;

   // this is not used, but it's needed in call to prtbp_del_2d
   double ti;   // integration time

   // Compute the 2D poincare map $P^(pt)$
   if(prtbp_del_2d(mu,sec,H,n,pt,&ti))
   {
      fprintf(stderr, "p12: error computing 2D poincare map\n");
      exit(EXIT_FAILURE);
   }
   return(pt[0]);
}

// Second component of Poincare map as a function of the first variable (g)
double p21(double g, void *params)
{
   double mu=((struct dparams *)params)->mu;
   section_t sec=((struct dparams *)params)->sec;
   double H=((struct dparams *)params)->H;
   int n=((struct dparams *)params)->cuts;
   double pt[2];

   pt[0]=g;
   pt[1]=((struct dparams *)params)->arg;	// G

   // this is not used, but it's needed in call to prtbp_del_2d
   double ti;   // integration time

   // Compute the 2D poincare map $P^(pt)$
   if(prtbp_del_2d(mu,sec,H,n,pt,&ti))
   {
      fprintf(stderr, "p21: error computing 2D poincare map\n");
      exit(EXIT_FAILURE);
   }
   return(pt[1]);
}

// Second component of Poincare map as a function of the second variable (G)
double p22(double G, void *params)
{
   double mu=((struct dparams *)params)->mu;
   section_t sec=((struct dparams *)params)->sec;
   double H=((struct dparams *)params)->H;
   int n=((struct dparams *)params)->cuts;
   double pt[2];

   pt[0]=((struct dparams *)params)->arg;	// g
   pt[1]=G;

   // this is not used, but it's needed in call to prtbp_del_2d
   double ti;   // integration time

   // Compute the 2D poincare map $P^(pt)$
   if(prtbp_del_2d(mu,sec,H,n,pt,&ti))
   {
      fprintf(stderr, "p22: error computing 2D poincare map\n");
      exit(EXIT_FAILURE);
   }
   return(pt[1]);
}

// First component of Poincare map as a function of the first variable (g)

double p11_inv(double g, void *params)
{
   double mu=((struct dparams *)params)->mu;
   section_t sec=((struct dparams *)params)->sec;
   double H=((struct dparams *)params)->H;
   int n=((struct dparams *)params)->cuts;
   double pt[2];

   pt[0]=g;
   pt[1]=((struct dparams *)params)->arg;	// G

   // this is not used, but it's needed in call to prtbp_del_2d
   double ti;   // integration time

   // Compute the 2D poincare map $P^(pt)$
   if(prtbp_del_2d_inv(mu,sec,H,n,pt,&ti))
   {
      fprintf(stderr, "p11_inv: error computing 2D poincare map\n");
      exit(EXIT_FAILURE);
   }
   return(pt[0]);
}

// First component of Poincare map as a function of the second variable (G)
double p12_inv(double G, void *params)
{
   double mu=((struct dparams *)params)->mu;
   section_t sec=((struct dparams *)params)->sec;
   double H=((struct dparams *)params)->H;
   int n=((struct dparams *)params)->cuts;
   double pt[2];

   pt[0]=((struct dparams *)params)->arg;	// g
   pt[1]=G;

   // this is not used, but it's needed in call to prtbp_del_2d
   double ti;   // integration time

   // Compute the 2D poincare map $P^(pt)$
   if(prtbp_del_2d_inv(mu,sec,H,n,pt,&ti))
   {
      fprintf(stderr, "p12_inv: error computing 2D poincare map\n");
      exit(EXIT_FAILURE);
   }
   return(pt[0]);
}

// Second component of Poincare map as a function of the first variable (g)
double p21_inv(double g, void *params)
{
   double mu=((struct dparams *)params)->mu;
   section_t sec=((struct dparams *)params)->sec;
   double H=((struct dparams *)params)->H;
   int n=((struct dparams *)params)->cuts;
   double pt[2];

   pt[0]=g;
   pt[1]=((struct dparams *)params)->arg;	// G

   // this is not used, but it's needed in call to prtbp_del_2d
   double ti;   // integration time

   // Compute the 2D poincare map $P^(pt)$
   if(prtbp_del_2d_inv(mu,sec,H,n,pt,&ti))
   {
      fprintf(stderr, "p21_inv: error computing 2D poincare map\n");
      exit(EXIT_FAILURE);
   }
   return(pt[1]);
}

// Second component of Poincare map as a function of the second variable (G)
double p22_inv(double G, void *params)
{
   double mu=((struct dparams *)params)->mu;
   section_t sec=((struct dparams *)params)->sec;
   double H=((struct dparams *)params)->H;
   int n=((struct dparams *)params)->cuts;
   double pt[2];

   pt[0]=((struct dparams *)params)->arg;	// g
   pt[1]=G;

   // this is not used, but it's needed in call to prtbp_del_2d
   double ti;   // integration time

   // Compute the 2D poincare map $P^(pt)$
   if(prtbp_del_2d_inv(mu,sec,H,n,pt,&ti))
   {
      fprintf(stderr, "p22_inv: error computing 2D poincare map\n");
      exit(EXIT_FAILURE);
   }
   return(pt[1]);
}

// name OF FUNCTION: dprtbp_del_2d
// CREDIT: 
//
// PURPOSE
// =======
// Consider the RTBP. Let "sec" be the Poincare section.
// Fix the Hamiltonian to a given value "H". 
// Using this energy condition, we can work with only two variables, $(g,G)$.
// The third variable $l$ is determined by Poincare section, and the
// fourth variable $L$ can be obtained from the energy condition.
// This function computes the derivative $DP(g,G)$ of the 2D Poincare map.
// The 2D derivative is computed NUMERICALLY from the 2D map using central
// differences.
//
// PARAMETERS
// ==========
// mu
//    mass parameter for the RTBP
//
// sec
//    type of Poincare section = {SEC1,SEC2}
//
// H
//    energy value
//
// cuts
//    number of iterates of the Poincare map (cuts=n: n cuts with the
//    Poincare section).
// 
// x
//    Argument of the 2d derivative, 2 coordinates: x=(g,G). 
//
// dp2d
//    On exit, dp2d holds the (2-by-2) derivative $DP$ of the 2D Poincare
//    map.
// 
// RETURN VALUE
// ============
// Returns a non-zero error code to indicate an error and 0 to indicate
// success.
//
// NOTES
// =====
//
// CALLS TO: p11,p12,p21,p22

int dprtbp_del_2d(double mu, section_t sec, double H, int cuts, double x[2], double dp2d[4])
{
   double abserr[4];		// error of numerical derivative
   struct dparams params;

   gsl_function P11, P12, P21, P22;

   // Compute derivative using numerical differentiation.
   params.mu = mu;
   params.sec = sec;
   params.H = H;
   params.cuts = cuts;

   params.arg = x[1];
   P11.function = &p11;
   P11.params = &params;
   gsl_deriv_central(&P11, x[0], STEP_SIZE, dp2d, abserr);

   params.arg = x[0];
   P12.function = &p12;
   P12.params = &params;
   gsl_deriv_central(&P12, x[1], STEP_SIZE, dp2d+1, abserr+1);

   params.arg = x[1];
   P21.function = &p21;
   P21.params = &params;
   gsl_deriv_central(&P21, x[0], STEP_SIZE, dp2d+2, abserr+2);

   params.arg = x[0];
   P22.function = &p22;
   P22.params = &params;
   gsl_deriv_central(&P22, x[1], STEP_SIZE, dp2d+3, abserr+3);

   fprintf(stderr, "Estimate of absolute error in numerical derivative:\n");
   fprintf(stderr, "% .10le % .10le\n", abserr[0], abserr[1]);
   fprintf(stderr, "% .10le % .10le\n", abserr[2], abserr[3]);
   fprintf(stderr, "\n");

   return(0);
}

// name OF FUNCTION: dprtbp_del_2d_inv
// CREDIT: 
//
// PURPOSE
// =======
// Consider the RTBP. Let S be the Poincare section {l=0}.
// Fix the Hamiltonian to a given value "H". 
// Using this energy condition, we can work with only two variables, $(g,G)$.
// The third variable $l$ is 0 since we look at the Poincare section, and the
// fourth variable $L$ can be obtained from the energy condition.
// This function computes the derivative $DP(g,G)$ of the 2D Poincare map.
// The 2D derivative is computed NUMERICALLY from the 2D map using central
// differences.
//
// PARAMETERS
// ==========
// mu
//    mass parameter for the RTBP
//
// H
//    energy value
//
// cuts
//    number of iterates of the Poincare map (cuts=n: n cuts with the
//    Poincare section).
// 
// x
//    Argument of the 2d derivative, 2 coordinates: x=(g,G). 
//
// dp2d
//    On exit, dp2d holds the (2-by-2) derivative $DP$ of the 2D Poincare
//    map.
// 
// RETURN VALUE
// ============
// Returns a non-zero error code to indicate an error and 0 to indicate
// success.
//
// NOTES
// =====
//
// CALLS TO: p11_inv,p12_inv,p21_inv,p22_inv

int dprtbp_del_2d_inv(double mu, double H, int cuts, double x[2], double dp2d[4])
{
   double abserr[4];		// error of numerical derivative
   struct dparams params;

   gsl_function P11, P12, P21, P22;

   // Compute derivative using numerical differentiation.
   params.mu = mu;
   params.H = H;
   params.cuts = cuts;

   params.arg = x[1];
   P11.function = &p11_inv;
   P11.params = &params;
   gsl_deriv_central(&P11, x[0], STEP_SIZE, dp2d, abserr);

   params.arg = x[0];
   P12.function = &p12_inv;
   P12.params = &params;
   gsl_deriv_central(&P12, x[1], STEP_SIZE, dp2d+1, abserr+1);

   params.arg = x[1];
   P21.function = &p21_inv;
   P21.params = &params;
   gsl_deriv_central(&P21, x[0], STEP_SIZE, dp2d+2, abserr+2);

   params.arg = x[0];
   P22.function = &p22_inv;
   P22.params = &params;
   gsl_deriv_central(&P22, x[1], STEP_SIZE, dp2d+3, abserr+3);

   fprintf(stderr, "Estimate of absolute error in numerical derivative:\n");
   fprintf(stderr, "% .10le % .10le\n", abserr[0], abserr[1]);
   fprintf(stderr, "% .10le % .10le\n", abserr[2], abserr[3]);
   fprintf(stderr, "\n");

   return(0);
}

