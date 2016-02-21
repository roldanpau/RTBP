// ===============================================================
// Derivative of Poincare map of the Restricted Three Body Problem
// ===============================================================
// FILE:          $RCSfile: dprtbp.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2012-11-26 09:48:57 $
//
// PURPOSE:
//
// NOTES:
//
// OVERALL METHOD:
//
// FUNCTIONS
// =========
//
// dprtbp
//    Compute the derivative of Poincare map $DP(x)$ of the RTBP.
//
// dprtbp_inv
//    Compute the derivative of inverse Poincare map $DP^{-1}(x)$ of the RTBP.

#include <stdio.h>	// fprintf
#include <prtbp.h>	// prtbp
#include <frtbp.h>	// dfrtbp
#include <rtbp.h>	// DIM

// name OF FUNCTION: dprtbp
// CREDIT: 
//
// DESCRIPTION
// ===========
// Compute the derivative of the Poincare map $DP^n(x)$ of the RTBP.
// Let $\phi(t,x)$ be the flow of the RTBP.
// The derivative of the Poincare map P in x is:
//    DP(x) = D\phi(T,x), 
// where T is the time to reach the poincare section.
// The derivative corresponds to the solution of the following first
// variational equation at T:
//    (D\phi(t,x))' = D_x F(t,x) D\phi(t,x),
// with initial condition
//    D\phi(t,x)=I.
// Integrating these first order variational equations from x along the
// trajectory which returns to the section, we obtain the derivative DP(x).
// 
// PARAMETERS
// ==========
// mu
//    mass parameter for the RTBP
// sec
//    Poincare section
// cuts
//    number of iterates of the Poincare map (cuts=n: n cuts with the
//    Poincare section).
// x
//    Argument of the derivative, 4 coordinates: (X, Y, P_X, P_Y).
// dp
//    On return of the this function, it holds the derivative $DP(x)$.
// 
// RETURN VALUE
// ============
// Returns a non-zero error code to indicate an error and 0 to indicate
// success.
// If an integration error is encountered, the function returns a non-zero
// value.
//
// NOTES
// =====
// On output, the point x is unmodified.
//
// CALLS TO: prtbp, dfrtbp

int dprtbp(double mu, section_t sec, int cuts, double x[DIM], double dp[DIMV])
{
   double t;	/* time to reach poincare section */
   int status,i;
   double x0[DIM];

   // Initial condition
   for(i=0; i<DIM; i++)
      x0[i]=x[i];

   // Compute time to reach poincare section
   status=prtbp(mu,sec,cuts,x0,&t);
   if(status)
   {
      fprintf(stderr, \
            "drtbp: error computing time to reach poincare section\n");
      return(1);
   }
   // Reset initial condition
   for(i=0; i<DIM; i++)
      x0[i]=x[i];

   // Integrate variational equations
   status=dfrtbp(mu,t,x0,dp);
   if(status)
   {
      fprintf(stderr, "drtbp: error integrating variational equations\n");
      return(1);
   }
   return(0);
}

// name OF FUNCTION: dprtbp_inv
// CREDIT: 
//
// DESCRIPTION
// ===========
// Compute the derivative of the inverse Poincare map $DP^{-n}(x)$ of the RTBP.
// Let $\phi(t,x)$ be the flow of the RTBP.
// The derivative of the inverse Poincare map P in x is:
//    DP^{-1}(x) = D\phi(T,x), 
// where T<0 is the time to reach the poincare section.
// The derivative corresponds to the solution of the following first
// variational equation at T<0:
//    (D\phi(t,x))' = D_x F(t,x) D\phi(t,x),
// with initial condition
//    D\phi(t,x)=I.
// Integrating these first order variational equations backwards from x along
// the trajectory which returns to the section, we obtain the derivative
// DP^{-1}(x).
// 
// PARAMETERS
// ==========
// mu
//    mass parameter for the RTBP
// sec
//    Poincare section
// cuts
//    number of iterates of the Poincare map (cuts=n: n cuts with the
//    Poincare section).
// x
//    Argument of the derivative, 4 coordinates: (X, Y, P_X, P_Y).
// dp
//    On return of the this function, it holds the derivative $DP^{-1}(x)$.
// 
// RETURN VALUE
// ============
// Returns a non-zero error code to indicate an error and 0 to indicate
// success.
// If an integration error is encountered, the function returns a non-zero
// value.
//
// NOTES
// =====
// On output, the point x is unmodified.
//
// CALLS TO: prtbp_inv, dfrtbp

int dprtbp_inv(double mu, section_t sec, int cuts, double x[DIM], double dp[DIMV])
{
   double t;	/* time to reach poincare section */
   int status,i;
   double x0[DIM];

   // Initial condition
   for(i=0; i<DIM; i++)
      x0[i]=x[i];

   // Compute time to reach poincare section. This time is negative.
   status=prtbp_inv(mu,sec,cuts,x0,&t);
   if(status)
   {
      fprintf(stderr, \
            "drtbp: error computing time to reach poincare section\n");
      return(1);
   }
   // Reset initial condition
   for(i=0; i<DIM; i++)
      x0[i]=x[i];

   // Integrate variational equations in negative time.
   status=dfrtbp(mu,t,x0,dp);
   if(status)
   {
      fprintf(stderr, "drtbp: error integrating variational equations\n");
      return(1);
   }
   return(0);
}
