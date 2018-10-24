/*! \file
    \brief Outer Map of the Circular Problem for Stochastic paper

    Compute the integrals $\omega_+^j$ and $\omega_-^j$ related to the outer
    map of the circular problem.
*/

#include <stdio.h>	// perror
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE
#include <string.h> // strcmp
#include <assert.h>

#include <gsl/gsl_integration.h>	// gsl_integration_qags
#include <utils_module.h>           // dblcpy
#include <section.h>
#include <rtbpdel.h>            	// f0_stoch
#include <frtbpred.h>
#include <prtbp_del_car.h>

// We request a absolute error of 0 and a relative error $10^{-13}$.

// NOTE: we can't reach accuracy of 10^{-13} when computing
// omega_pos_f, so we lower it to 10^{-9}
const double INTEGRATION_EPSABS = 0.0;
const double INTEGRATION_EPSREL = 1.e-8;

/// Parameters to the \ref integrand_omega_pm function.
struct iparams_omega_pm
{
   double mu;
   double x[DIM];
};

double integrand_omega_pm(double s, void *params)
{
   double mu;
   double x[DIMRED];

   // auxiliary variables
   int i,status;
   double x2[DIM];
   double res;      // res=f0(\gamma(s))

   mu = ((struct iparams_omega_pm *)params)->mu;

   x[0] = ((struct iparams_omega_pm *)params)->x[0];    // l
   x[1] = ((struct iparams_omega_pm *)params)->x[1];    // L
   x[2] = ((struct iparams_omega_pm *)params)->x[2];    // g
   x[3] = ((struct iparams_omega_pm *)params)->x[3];    // G

   x[4] = 0;    // t
   x[5] = 0;    // I (not used).

   // Compute x = \lambda(s)
   status = frtbp_red_g(mu,s,x);
   if(status)
   {
      fprintf(stderr, "integrand_omega_pm: error integrating trajectory");
      exit(EXIT_FAILURE);
   }

   // Compute f0(x)
   for(i=0;i<DIM;i++) x2[i]=x[i];
   status = f0_stoch(x2,&res,&mu);
   if(status)
   {
      fprintf(stderr, "integrand_omega_pm: error computing f0");
      exit(EXIT_FAILURE);
   }
   return res;
}


// NOTE: Instead of P^{-(N-i)}(z^s), we could have used
// frtbp_red(-2(N-i)\pi, z^s). They should give the same point.

int omega_pos_stoch(double mu, section_t sec, double x[DIM], double x_car[DIM],
		int N, double T0, double *omega) 
{
   double result, error;

   // auxiliary variables
   int i,j;
   double xi[DIM];		    /* point \xi=P^{-(N-i)}(z^s) in Delaunay */
   double xi_car[DIM];		/* point \xi=P^{-(N-i)}(z^s) in Cartesian */
   double t;

   assert(N>0);

   gsl_integration_workspace * w
               = gsl_integration_workspace_alloc (1000);

   struct iparams_omega_pm params;

   gsl_function F;
   F.function = &integrand_omega_pm;
   F.params = &params;

   params.mu = mu;
   
   *omega = 0.0;
   for(i=N; i>=1; i--)
   {
      // $\gamma_i(s)$ is the homoclinic trajectory that starts at the
      // point \xi = P^{-(N-i)}(z^s) = P^{i}(z). 
       dblcpy(xi,x,DIM);
       dblcpy(xi_car,x_car,DIM);

      if(prtbp_del_car_inv(mu,sec,(N-i),xi,xi_car,&t))
      {
         fprintf(stderr, 
                 "omega_pos_stoch: error computing point P^{-%d}(z^s)\n", 
                 N-i);
         return(1);
      }

      dblcpy(params.x, xi, DIM);

      // Integrate integrand function from -2\pi to 0. 
      gsl_integration_qags (&F, -2*M_PI, 0, INTEGRATION_EPSABS,
			  INTEGRATION_EPSREL, 1000, w, &result, &error); 
      fprintf (stderr, "estimated error = % .3le\n", error);

      // \omega_+
      (*omega) = (*omega) +(result + T0);
   }

   gsl_integration_workspace_free (w);
   return 0;
}

// NOTE: Instead of P^{N-i}(z^u), we could have used
// frtbp_red(2(N-i)\pi, z^u). They should give the same point.

int omega_neg_stoch(double mu, section_t sec, double x[DIM], double x_car[DIM], 
		int N, double T0, double *omega)
{
   double result, error;

   // auxiliary variables
   int i,j;
   double xi[DIM];		    /* point \xi=P^{N-i}(z^u) */
   double xi_car[DIM];		/* point \xi=P^{N-i}(z^u) */
   double t;

   assert(N>0);

   gsl_integration_workspace * w
               = gsl_integration_workspace_alloc (1000);

   struct iparams_omega_pm params;

   gsl_function F;
   F.function = &integrand_omega_pm;
   F.params = &params;

   params.mu = mu;
   
   *omega = 0.0;
   for(i=N; i>=1; i--)
   {
      // $\gamma_i(s)$ is the homoclinic trajectory that starts at the
      // point \xi = P^{N-i}(z^u) = P^{-i}(z). 
       dblcpy(xi,x,DIM);
       dblcpy(xi_car,x_car,DIM);

      if(prtbp_del_car(mu,sec,(N-i),xi,xi_car,&t))
      {
         fprintf(stderr, "omega_neg_stoch: error computing point P^{%d}(z^u)\n",
               (N-i));
         return(1);
      }

      dblcpy(params.x, xi, DIM);

      // Integrate integrand function from 2\pi to 0. 
      gsl_integration_qags (&F, 2*M_PI, 0, INTEGRATION_EPSABS,
			  INTEGRATION_EPSREL, 1000, w, &result, &error); 
      fprintf (stderr, "estimated error = % .3le\n", error);

      // \omega_-
      (*omega) = (*omega) +(result - T0);
   }

   gsl_integration_workspace_free (w);
   return 0;
}
