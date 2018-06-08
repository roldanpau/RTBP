/*! \file 

  \brief Inner Map of the Elliptic Problem for Stochastic paper

  Compute the integrals $B_{in}^j$ for $j=1,2$ related to the inner map
  of the elliptic problem.

  This program computes the shift of the inner map of the elliptic problem,
  given by the integral
  
  \[ A_1(I,t) = A_1^+ e^{it} + A_1^- e^{-it}. \]  ???
 
  \[ A_{in}^j(I) = \mu\int_0^{2\pi} f(\lambda(s)) e^{it(s)} ds, \]
  where
  \[ f(l,L,g,G) = \frac{\Delta H_{ell}^{1,+}(l,L,g,G)}
        {-1+\mu\partial_G \Delta H_{circ}}(l,L,g,G). \]
  and $\lambda(s)$ is the trajectory of the reduced circular problem with
  initial condition (l,L,0,G,t=0).
  
  The integrand is evaluated at an equispaced sequence of points 
  t\in(0,2\pi). 
*/

// FUNCTIONS
// =========
//
// re_f_integrand_stoch
// im_f_integrand_stoch
// --------------
// The function
// \[ f(l,L,g,G) = \frac{\Delta H_{ell}^{1,+}(l,L,g,G)}
//       {-1+\mu\partial_G \Delta H_{circ}}(l,L,g,G). \]
//
// re_integrand_inner_ell_stoch
// im_integrand_inner_ell_stoch
// ----------------------
// Given the unique point (l,L,g=0,G) that belongs to the intersection
// between the periodic orbit and with the section $g=0$,
// and an integration time 's' in the periodic trajectory, 
// this function computes the integrand $f(\gamma(s)) e^{it(s)}$.
//
// re_inner_ell_stoch
// im_inner_ell_stoch
// ------------
// Given an energy level $H$, compute the integral $A_{in}^f$.
// This is computed using numerical integration.

#include <stdio.h>	// perror
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE

#include <gsl/gsl_integration.h>	// gsl_integration_qags

#include <utils_module.h>	// dblcpy
#include <rtbpdel.h>		// dot_g, re_DHell, im_DHell
#include <frtbpred.h>

struct iparams_inner_ell_stoch
{
   double mu;
   double x[DIM];
};

// name OF FUNCTION: re_f_integrand_stoch
//
// PURPOSE
// =======
// Real part of the function
// \[ f(l,L,g,G) = \frac{\Delta H_{ell}^{1,+}(l,L,g,G)}
//       {-1+\mu\partial_G \Delta H_{circ}}(l,L,g,G). \]
//
// PARAMETERS
// ==========
// mu
//    mass parameter for the RTBP
// x
//    point, 4 coordinates: (l,L,g,G). 
// 
// RETURN VALUE
// ============
// Returns the REAL part of the function $f$.
//
// NOTES
// =====
//
// CALLS TO: dot_g, re_DHell

double re_f_integrand_stoch(double mu, double x2[DIM])
{
   // auxiliary variables
   int status;
   double y[DIM];

   double den;		// denominator of f
   double re_num;	// numerator of f (real part)

   // Compute the denominator $-1+\mu\partial_G \Delta H_{circ}}$.
   // This is just the $\dot g$ component in the nonreduced vector field.
   status = dot_g(x2,&den,&mu);
   if(status)
   {
      fprintf(stderr, "re_f_integrand_stoch: error computing denominator");
      exit(EXIT_FAILURE);
   }

   // Compute the real part of the numerator 
   re_num = re_DHell(x2,&mu);

   return re_num/den;
}

// name OF FUNCTION: im_f_integrand_stoch
//
// PURPOSE
// =======
// Imaginary part of the function
// \[ f(l,L,g,G) = \frac{\Delta H_{ell}^{1,+}(l,L,g,G)}
//       {-1+\mu\partial_G \Delta H_{circ}}(l,L,g,G). \]
//
// PARAMETERS
// ==========
// mu
//    mass parameter for the RTBP
// x
//    point, 4 coordinates: (l,L,g,G). 
// 
// RETURN VALUE
// ============
// Returns the REAL part of the function $f$.
//
// NOTES
// =====
//
// CALLS TO: dot_g, im_DHell

double im_f_integrand_stoch(double mu, double x2[DIM])
{
   // auxiliary variables
   int status;
   double y[DIM];

   double den;		// denominator of f
   double im_num;	// numerator of f (imaginary part)

   // Compute the denominator $-1+\mu\partial_G \Delta H_{circ}}$.
   // This is just the $\dot g$ component in the nonreduced vector field.
   status = dot_g(x2,&den,&mu);
   if(status)
   {
      fprintf(stderr, "im_f_integrand_stoch: error computing denominator");
      exit(EXIT_FAILURE);
   }

   // Compute the imaginary part of the numerator 
   im_num = im_DHell(x2,&mu);

   return im_num/den;
}

// name OF FUNCTION: re_integrand_inner_ell_stoch
//
// PURPOSE
// =======
// Given the unique point x=(l,L,g=0,G) that belongs to the intersection
// between the periodic orbit and with the section $g=0$,
// and an integration time 's' in the periodic trajectory, 
// this function computes the integrand $f(\gamma(s)) e^{it(s)}$.
//
// PARAMETERS
// ==========
// s
//    integration time in the periodic trajectory.
// mu
//    mass parameter for the RTBP
// x
//    periodic point, 4 coordinates: (l,L,g,G). 
// 
// RETURN VALUE
// ============
// Returns the REAL part of the integrand evaluated at the point $\lambda(s)$.
//
// NOTES
// =====
//
// CALLS TO: frtbp_red_g, re_f_integrand_stoch, im_f_integrand_stoch

double re_integrand_inner_ell_stoch(double s, void *params)
{
   double mu;
   double x[DIMRED];

   // auxiliary variables
   int status;
   double x2[DIM], y[DIM];

   double re_f, im_f;	// f(\lambda(s))
   double t;		// original time

   mu = ((struct iparams_inner_ell_stoch *)params)->mu;

   x[0] = ((struct iparams_inner_ell_stoch *)params)->x[0]; 	// l
   x[1] = ((struct iparams_inner_ell_stoch *)params)->x[1]; 	// L
   x[2] = ((struct iparams_inner_ell_stoch *)params)->x[2]; 	// g
   x[3] = ((struct iparams_inner_ell_stoch *)params)->x[3]; 	// G

   x[4] = 0;	// t
   x[5] = 0;	// I (not used).

   // Compute x = \lambda(s)
   status = frtbp_red_g(mu,s,x);
   if(status)
   {
      fprintf(stderr, "integrand: error integrating trajectory");
      exit(EXIT_FAILURE);
   }
   t = x[4];

   // x2=\lambda(s). Drop the last two components of x.
   dblcpy(x2,x,DIM);

   // Compute the function f(\lambda(s))
   re_f = re_f_integrand_stoch(mu,x2);
   im_f = im_f_integrand_stoch(mu,x2);

   // PRG (6/6/18): I believe this was computed wrong in the last paper!!!
   //return re_f*cos(t) + im_f*sin(t);
   return -(re_f*sin(t) + im_f*cos(t));
}

// name OF FUNCTION: im_integrand_inner_ell_stoch
//
// PURPOSE
// =======
// Given the unique point x=(l,L,g=0,G) that belongs to the intersection
// between the periodic orbit and with the section $g=0$,
// and an integration time 's' in the periodic trajectory, 
// this function computes the integrand $f(\lambda(s)) e^{it(s)}$.
//
// PARAMETERS
// ==========
// s
//    integration time in the periodic trajectory.
// mu
//    mass parameter for the RTBP
// x
//    periodic point, 4 coordinates: (l,L,g,G). 
// 
// RETURN VALUE
// ============
// Returns the IMAGINARY part of the integrand evaluated at the point
// $\lambda(s)$.
//
// NOTES
// =====
//
// CALLS TO: frtbp_red_g, re_f_integrand_stoch, im_f_integrand_stoch

double im_integrand_inner_ell_stoch(double s, void *params)
{
   double mu;
   double x[DIMRED];

   // auxiliary variables
   int status;
   double x2[DIM], y[DIM];

   double re_f, im_f;	// f(\lambda(s))
   double t;		// original time

   mu = ((struct iparams_inner_ell_stoch *)params)->mu;

   x[0] = ((struct iparams_inner_ell_stoch *)params)->x[0]; 	// l
   x[1] = ((struct iparams_inner_ell_stoch *)params)->x[1]; 	// L
   x[2] = ((struct iparams_inner_ell_stoch *)params)->x[2]; 	// g
   x[3] = ((struct iparams_inner_ell_stoch *)params)->x[3]; 	// G

   x[4] = 0;	// t
   x[5] = 0;	// I (not used).

   // Compute x = \lambda(s)
   status = frtbp_red_g(mu,s,x);
   if(status)
   {
      fprintf(stderr, "integrand: error integrating trajectory");
      exit(EXIT_FAILURE);
   }
   t = x[4];

   // x2=\lambda(s). Drop the last two components of x.
   dblcpy(x2,x,DIM);

   // Compute the function f(\lambda(s))
   re_f = re_f_integrand_stoch(mu,x2);
   im_f = im_f_integrand_stoch(mu,x2);

   // PRG (6/6/18): I believe this was computed wrong in the last paper!!!
   //return re_f*sin(t) - im_f*cos(t);
   return re_f*cos(t) - im_f*sin(t);
}

// name OF FUNCTION: re_inner_ell_stoch
// CREDIT: 
//
// PURPOSE
// =======
// The inner map for the elliptic problem is given by
// 
//    \[ B_{in}^j(I) = i\mu \frac{1-e^{i\omega}}{1-e^{i2\pi\nu}} 
//       \int_0^{T} \frac{\Delta H_{ell}^{1,+}}
//       {-1+\mu\partial_G \Delta H_{circ}} e^{it(s)} ds. \]
//
// Given an energy level $H$, this function computes the real part of
// \[ A_{in}(I) := i\mu \int_0^{T} \frac{\Delta H_{ell}^{1,+}} 
//       {-1+\mu\partial_G \Delta H_{circ}} e^{it(s)} ds. \]
//
// PARAMETERS
// ==========
// mu
//    mass parameter for the RTBP
// T
//    limit of integration. This will be 2\pi for B_{in}^j.
// x
//    x=(l,L,g=0,G), periodic point of period 1, on the section g=0.
// re_A
//    On return of this function, re_A contains the real part of $A_{in}(I)$.
// 
// RETURN VALUE
// ============
// Returns a non-zero error code to indicate an error and 0 to indicate
// success.
//
// CALLS TO: re_integrand_inner_ell_stoch

int re_inner_ell_stoch(double mu, double T, double x[DIM], double *re_A)
{
   double result, error;
   gsl_integration_workspace * w
               = gsl_integration_workspace_alloc (1000);

   struct iparams_inner_ell_stoch params;

   params.mu = mu;
   
   params.x[0] = x[0];
   params.x[1] = x[1];
   params.x[2] = x[2];
   params.x[3] = x[3];

   gsl_function F;
   F.function = &re_integrand_inner_ell_stoch;
   F.params = &params;

   // Integrate integrand function from 0 to T. 
   // We request a absolute error of 0 and a relative error $10^{-9}$.
   // NOTE: this relative error is the same as the one used for outer_circ.
   gsl_integration_qags (&F, 0, T, 0, 1.e-9, 1000, w, &result, &error);
   fprintf (stderr, "estimated error = % .3le\n", error);

   gsl_integration_workspace_free (w);

   *re_A = mu*result;		// real(A^+)
   return 0;
}

// name OF FUNCTION: im_inner_ell_stoch
// CREDIT: 
//
// PURPOSE
// =======
// Given an energy level $H$, this function computes the imaginary part of
// \[ A_{in}(I) := i\mu \int_0^{T} \frac{\Delta H_{ell}^{1,+}} 
//       {-1+\mu\partial_G \Delta H_{circ}} e^{it(s)} ds. \]

int im_inner_ell_stoch(double mu, double T, double x[DIM], double *im_A)
{
   double result, error;
   gsl_integration_workspace * w
               = gsl_integration_workspace_alloc (1000);

   struct iparams_inner_ell_stoch params;

   params.mu = mu;
   
   params.x[0] = x[0];
   params.x[1] = x[1];
   params.x[2] = x[2];
   params.x[3] = x[3];

   gsl_function F;
   F.function = &im_integrand_inner_ell_stoch;
   F.params = &params;

   // Integrate integrand function from 0 to T. 
   // We request a absolute error of 0 and a relative error $10^{-9}$.
   // NOTE: this relative error is the same as the one used for outer_circ.
   gsl_integration_qags (&F, 0, T, 0, 1.e-9, 1000, w, &result, &error);
   fprintf (stderr, "estimated error = % .3le\n", error);

   gsl_integration_workspace_free (w);

   *im_A = mu*result;		// imaginary(A^+)
   return 0;
}
