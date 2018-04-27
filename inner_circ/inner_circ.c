/*! \file
    \brief Inner Map of the Circular Problem.

    $Author: roldan $
    $Date: 2013-03-26 22:20:08 $
*/

#include <stdio.h>	// perror
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE

#include <gsl/gsl_integration.h>	// gsl_integration_qags
#include <rtbpdel.h>			// dot_g, f0
#include <frtbpred.h>

#include "inner_circ.h" 	// iparams_omega_in

// FUNCTIONS
// =========
//
// integrand_inner_circ
// --------------------
// Given a periodic trajectory \gamma(s) of the circular RTBP, 
// and an integration time 's' in the periodic trajectory, 
// this function computes the integrand,
//    \frac{1}{-1+\mu\partial_G \Delta H_{circ}(\gamma(s))}.
//
// NOTE: This function is not used anymore, since $T_0$ (inner map of circular
// problem) is actually computed from the period of the periodic orbit:
//
// \[ 2\pi + \mu T_0 = period of p.o. \]
//

// Parameters to \ref inner_circ function.
struct iparams_inner_circ
{
   double mu;
   double x[DIM];
};

double integrand_omega_in(double s, void *params)
{
   double mu;
   double x[DIMRED];

   // auxiliary variables
   int i,status;
   double x2[DIM];
   double res;		// res=f0(\gamma(s))

   mu = ((struct iparams_omega_in *)params)->mu;

   x[0] = ((struct iparams_omega_in *)params)->x[0]; 	// l
   x[1] = ((struct iparams_omega_in *)params)->x[1]; 	// L
   x[2] = ((struct iparams_omega_in *)params)->x[2]; 	// g
   x[3] = ((struct iparams_omega_in *)params)->x[3]; 	// G

   x[4] = 0;	// t
   x[5] = 0;	// I (not used).

   // Compute x = \lambda(s)
   status = frtbp_red_l(mu,s,x);
   if(status)
   {
      fprintf(stderr, "integrand_omega_in: error integrating trajectory");
      exit(EXIT_FAILURE);
   }

   // Compute f0(x)
   for(i=0;i<DIM;i++) x2[i]=x[i];
   status = f0(x2,&res,&mu);
   if(status)
   {
      fprintf(stderr, "integrand_omega_in: error computing f0");
      exit(EXIT_FAILURE);
   }
   return res;
}

int omega_in_f(double mu, double x[DIM], double *omega)
{
   double result, error;
   gsl_integration_workspace * w
               = gsl_integration_workspace_alloc (1000);

   struct iparams_omega_in params;

   params.mu = mu;
   
   params.x[0] = x[0];
   params.x[1] = x[1];
   params.x[2] = x[2];
   params.x[3] = x[3];

   gsl_function F;
   F.function = &integrand_omega_in;
   F.params = &params;

   // Integrate integrand function from 0 to 4\pi. 
   // We request a absolute error of 0 and a relative error $10^{-13}$.
   // We request a absolute error of 0 and a relative error $10^{-10}$, since
   // 10^{-13} seems to be too much...
   gsl_integration_qags (&F, 0, 4*M_PI, 0, 1.e-10, 1000, w, &result, &error);
   fprintf (stderr, "estimated error = % .3le\n", error);

   gsl_integration_workspace_free (w);

   *omega = result;		// \omega_in^f
   return 0;
}

int omega_in_b(double mu, double x[DIM], double *omega)
{
   double result, error;
   gsl_integration_workspace * w
               = gsl_integration_workspace_alloc (1000);

   struct iparams_omega_in params;

   params.mu = mu;
   
   params.x[0] = x[0];
   params.x[1] = x[1];
   params.x[2] = x[2];
   params.x[3] = x[3];

   gsl_function F;
   F.function = &integrand_omega_in;
   F.params = &params;

   // Integrate integrand function from 0 to 2\pi. 
   // We request a absolute error of 0 and a relative error $10^{-13}$.
   // We request a absolute error of 0 and a relative error $10^{-10}$, since
   // 10^{-13} seems to be too much...
   gsl_integration_qags (&F, 0, 2*M_PI, 0, 1.e-10, 1000, w, &result, &error);
   fprintf (stderr, "estimated error = % .3le\n", error);

   gsl_integration_workspace_free (w);

   *omega = result;		// \omega_in^b
   return 0;
}

// name OF FUNCTION: integrand_inner_circ
//
// PURPOSE
// =======
// Let 
// 
// \[ I = 
//       \int_0^{6\pi} 
//          \frac{1}{L^{-3}+\mu\partial_L \Delta H_{circ}(\gamma(s))} ds, \]
//
// where $\gamma(s)$ is the periodic trajectory in the level of energy H.
// 
// Given a periodic trajectory \gamma(s) of the reduced circular RTBP, 
// and an integration time 's' in the periodic trajectory, 
// this function computes the integrand
//
//    \frac{1}{L^{-3}+\mu\partial_L \Delta H_{circ}(\gamma(s))}.
//
// NOTE: This function is not used anymore, since $T_0$ (inner map of circular
// problem) is actually computed from the period of the periodic orbit:
//
// \[ 2\pi + \mu T_0 = period of p.o. \]
//
// PARAMETERS
// ==========
// s
//    integration time in the periodic trajectory.
// mu
//    mass parameter for the RTBP
// x
//    periodic point, 4 coordinates:
//    (l,L,g,G). 
// 
// RETURN VALUE
// ============
// Returns the integrand $\frac{1}{L^{-3}+\mu\partial_L \Delta H_{circ}}$
// evaluated at the point $\gamma(s)$.
//
// NOTES
// =====
//
// CALLS TO: frtbp_red, dot_l

/**
  Integrand of inner circular problem.
  */
double integrand_inner_circ(double s, void *params)
{
   double mu;
   double x[DIMRED];

   // auxiliary variables
   int i,status;
   double x2[DIM];
   double dl;		// \dot l

   mu = ((struct iparams_inner_circ *)params)->mu;

   x[0] = ((struct iparams_inner_circ *)params)->x[0]; 	// l
   x[1] = ((struct iparams_inner_circ *)params)->x[1]; 	// L
   x[2] = ((struct iparams_inner_circ *)params)->x[2]; 	// g
   x[3] = ((struct iparams_inner_circ *)params)->x[3]; 	// G

   x[4] = 0;	// t
   x[5] = 0;	// I (not used).

   // Compute x = \gamma(s)
   status = frtbp_red_l(mu,s,x);
   if(status)
   {
      fprintf(stderr, "integrand: error integrating trajectory");
      exit(EXIT_FAILURE);
   }

   // Compute the denominator $L^{-3}+\mu\partial_L \Delta H_{circ}}$.
   // This is just the $\dot l$ component in the nonreduced vector field.
   for(i=0;i<DIM;i++) x2[i]=x[i];
   status = dot_l(x2,&dl,&mu);
   if(status)
   {
      fprintf(stderr, "integrand_inner_circ: error computing dot_l");
      exit(EXIT_FAILURE);
   }
   return 1.0/dl;
}

// NOTE: This function is not used anymore, since $T_0$ (inner map of circular
// problem) is actually computed from the period of the periodic orbit:
//
// \[ 2\pi + \mu T_0 = period of p.o. \]

// NOTE: What we denote $T_0$ here is sometimes refered as $\mu T_0$ in
// Marcel's notes.

int inner_circ(double mu, double x[DIM], double *T)
{
   double result, error;
   gsl_integration_workspace * w
               = gsl_integration_workspace_alloc (1000);

   struct iparams_inner_circ params;

   params.mu = mu;
   
   params.x[0] = x[0];
   params.x[1] = x[1];
   params.x[2] = x[2];
   params.x[3] = x[3];

   gsl_function F;
   F.function = &integrand_inner_circ;
   F.params = &params;

   // Integrate integrand function from 0 to 6\pi. 
   // We request a absolute error of 0 and a relative error $10^{-13}$.
   gsl_integration_qags (&F, 0, 6*M_PI, 0, 1.e-13, 1000, w, &result, &error);
   fprintf (stderr, "estimated error = % .3le\n", error);

   gsl_integration_workspace_free (w);

   // WARNING! this returns \mu T_0, not T_0!
   *T = result-2*M_PI;		// T_0
   return 0;
}
