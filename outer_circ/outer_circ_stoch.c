/*! \file
    \brief Outer Map of the Circular Problem for Stochastic paper: main prog.

    Compute the integrals $\omega_+^j$ and $\omega_-^j$ related to the outer
    map of the circular problem.


    $Author: roldan $
    $Date: 2013-03-26 22:32:02 $
*/

#include <stdio.h>	// perror
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE
#include <string.h> // strcmp
#include <assert.h>

#include <gsl/gsl_integration.h>	// gsl_integration_qags
#include <section.h>
#include <rtbpdel.h>            	// f0_stoch
#include <frtbpred.h>
#include <prtbp_del_car.h>

// We request a absolute error of 0 and a relative error $10^{-13}$.

// NOTE: we can't reach accuracy of 10^{-13} when computing
// omega_pos_f, so we lower it to 10^{-9}
const double INTEGRATION_EPSABS = 0.0;
const double INTEGRATION_EPSREL = 1.e-9;

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


/**
  Given an energy level \f$H\f$, compute \f$\omega_+^j(H)\f$.

  Given an energy level \f$H\f$, compute \f$\omega_+^j(H)\f$, defined as
 
  \f[ \omega_+^j(I) = 
        lim_{N\to \infty} 
           \int_0^{2\pi N} f0(\gamma^*(s)) ds + NT_0(I),
  \f]
 
  where \f$\gamma^*(s)\f$ is the homoclinic trajectory of the reduced flow that
  starts at the primary homoclinic point z.
 
  Equivalently, we compute it as
  
  \f[\omega_+^j(I) = 
        lim_{N\to +\infty} 
           \int_{-2\pi N}^0 f0(\gamma^*(s)) ds + NT_0(I),
  \f]
 
  where \f$\gamma^*(s)\f$ is the homoclinic trajectory that starts at the point 
  \f$ z^s = \mathcal{P}^{N}(z) \f$. Notice that z^s is in the STABLE manifold of
  the reduced flow.
 
  Equivalently, we split the integral from -2\pi N to 0 into N parts of size
  2\pi:
 
  \f[\omega_+^j(I) = 
        lim_{N\to +\infty} 
           \sum_{i=N,1} 
              (\int_{-2\pi}^{0} f0(\gamma_i(s)) ds + T_0(I)),
  \f]
 
  where \f$\gamma_i(s)\f$ is the homoclinic trajectory that starts at the point
  P^{-(N-i)}(z^s) = P^{i}(z). 
 
  This is computed using numerical integration.
 
  \param[in] mu 	mass parameter for the RTBP
  \param[in] sec    Poincare section: sec={SECg,SECg2}

  \param[in] x    [DIM]	x=(l,L,g,G),     point z^s, on the section g=0.
  \param[in] x_car[DIM]	x=(x,y,p_x,p_y), point z^s, on the section g=0.

  \param[in] N
     number of iterates of the poincare map (length of integration).
     N must be positive since we are computing \f$\omega_+^j\f$.

  \param[in] T0	shift of inner map

  \param[out] omega
     On return of this function, omega contains the value of the function
     \f$\omega_+^j\f$.
  
  \returns
  a non-zero error code to indicate an error and 0 to indicate success.
 */

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
      for(j=0;j<DIM;j++) xi[j]=x[j];
      for(j=0;j<DIM;j++) xi_car[j]=x_car[j];

      if(prtbp_del_car_inv(mu,sec,(N-i),xi,xi_car,&t))
      {
         fprintf(stderr, "omega_pos_stoch: error computing point P^{i}(z)\n");
         return(1);
      }

      params.x[0] = xi[0];
      params.x[1] = xi[1];
      params.x[2] = xi[2];
      params.x[3] = xi[3];

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

/**
  Given an energy level \f$H\f$, compute \f$\omega_-^j(H)\f$.

  Given an energy level \f$H\f$, compute \f$\omega_-^j(H)\f$, defined as
 
  \f[ \omega_-^j(I) = 
        lim_{N\to -\infty} 
           \left( \int_0^{2\pi N} f0(\gamma^*(s)) ds + NT_0(I) \right),
  \f]
 
  where \f$\gamma^*(s)\f$ is the homoclinic trajectory of the reduced flow that
  starts at the primary homoclinic point z.
 
  Equivalently, we compute it as
  
  \f[\omega_-^j(I) = 
        lim_{N\to \infty} 
           \left( \int_{2\pi N}^0 f0(\gamma^*(s)) ds - NT_0(I) \right),
  \f]
 
  where \f$\gamma^*(s)\f$ is the homoclinic trajectory that starts at the point 
  \f$ z^u = \mathcal{P}^{-N}(z) \f$. Notice that z^u is in the UNSTABLE
  manifold of the reduced flow.
 
  Equivalently, we split the integral from 2\pi N to 0 into N parts of size
  2\pi:
 
  \f[\omega_-^j(I) = 
        lim_{N\to \infty} 
           \left( \sum_{i=N,1} 
              (\int_{2\pi}^{0} f0(\gamma_i(s)) ds - T_0(I)) \right),
  \f]
 
  where \f$\gamma_i(s)\f$ is the homoclinic trajectory that starts at the point
  P^{N-i}(z^u) = P^{-i}(z). 
 
  This is computed using numerical integration.
 
  \param[in] mu 	mass parameter for the RTBP
  \param[in] sec    Poincare section: sec={SECg,SECg2}

  \param[in] x    [DIM]	x=(l,L,g,G),     point z^u, on the section g=0.
  \param[in] x_car[DIM]	x=(x,y,p_x,p_y), point z^u, on the section g=0.

  \param[in] N
     number of iterates of the poincare map (length of integration).
     N must be a POSITIVE integer.

  \param[in] T0	shift of inner map

  \param[out] omega
     On return of this function, omega contains the value of the function
     \f$\omega_-^j\f$.
  
  \returns
  a non-zero error code to indicate an error and 0 to indicate success.
 */

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
      for(j=0;j<DIM;j++) xi[j]=x[j];
      for(j=0;j<DIM;j++) xi_car[j]=x_car[j];

      if(prtbp_del_car(mu,sec,(N-i),xi,xi_car,&t))
      {
         fprintf(stderr, "omega_neg_stoch: error computing point P^{%d}(z^u)\n",
               (N-i));
         return(1);
      }

      params.x[0] = xi[0];
      params.x[1] = xi[1];
      params.x[2] = xi[2];
      params.x[3] = xi[3];

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

/**
  Outer Map of the Circular Problem: main prog.

  This program computes the integrals \f$ \omega_{+,-}^j \f$ related to the
  outer map of the circular problem.
  (See the new paper ``Stochastic Diffusion''.)

  It reads the following input from stdin:
  - mu
     mass parameter for the RTBP
  - sec
     Poincare section

  And a sequence of lines:
  - T
     period of periodic orbit.
  - zu
     preimage of homoclinic point z (in the original flow): P^{M}(z_u) = z.
     (Delaunay coordinates)
  - zu_car
     preimage of homoclinic point z (in the original flow): P^{M}(z_u) = z.
     (Cartesian coordinates)
  - M
     number of poincare iterates to reach z from z_u.

  For each input line, it outputs result to stdout:
  - omega_neg
 
 */
 
//  NOTE: Due to symmetry, we have the following relation: 
//     \omega_- = - \omega_+, 
//  so it is enough to compute one of them.
 
int main( )
{
   double mu;
   section_t sec;		/* Poincare section */

   double zu[DIM];	    /* preimage of primary homoclinic point */
   double zu_car[DIM];	/* preimage of primary homoclinic point (Cartesian) */

   double w_pos;	/* value of integral \omega_+^* */
   double w_neg;	/* value of integral \omega_-^* */
   double w_out;	/* value of integral \omega_{out}^* */

   double T, T0;		/* period,shift of inner map */

   /* Number of iterates to hit z from z_u: P^M(z_u)=z. 
      This will also be used as N, the number of iterates of poincare map
      along homoclinic orbit (length of integration) */
   int M;	

   // auxiliary vars
   char section_str[10];    // holds input string "SEC1", "SEC2" etc

   int i;
   double t;
   double w_neg_test;	/* value of integral \omega_-^* */

   // Input parameters from stdin.
   if(scanf("%le %s", &mu, section_str)<2)
   {
      perror("main: error reading input");
      exit(EXIT_FAILURE);
   }

   if (strcmp(section_str,"SEC1") == 0)
      sec = SEC1;
   else if (strcmp(section_str,"SEC2") == 0)
      sec = SEC2;
   else if (strcmp(section_str,"SECg") == 0)
      sec = SECg;
   else if (strcmp(section_str,"SECg2") == 0)
      sec = SECg2;
   else
   {
      perror("main: error reading section string");
      exit(EXIT_FAILURE);
   }

   // Input period T, zu, number of poincare iterates M, from stdin.
   while(scanf("%le %le %le %le %le %le %le %le %le %d", &T, 
               zu, zu+1, zu+2, zu+3, 
               zu_car, zu_car+1, zu_car+2, zu_car+3, 
               &M) == 10)
   {
      // TESTING...
      //T0 = (T-2*M_PI)/mu;
      T0 = T-2*M_PI;

      //omega_pos(mu, sec, zu, M, T0, &w_pos);

      // Compute $\omega_-^*$, integrating along $z(s) = \gamma^*(s)$.
      // Note: since the homoclinic point is at the symmetry axis, we have
      // \omega_-^* = -\omega_+^*.
      //w_neg = -w_pos;

      // TESTING
	  for(i=1; i<=M; i++) 
	  {
		  omega_neg_stoch(mu, sec, zu, zu_car, i, T0, &w_neg_test);
		  printf("%d %.15e \n", i, w_neg_test);

		  //omega_neg_stoch(mu, sec, zu, zu_car, M, T0, &w_neg);

		  //w_out = w_pos-w_neg;

		  // Output result to stdout.
		  //printf("%.15e\n", w_neg);
		  //printf("%.15e\n", w_pos);
		  fflush(NULL);
	  }
   }
   exit(EXIT_SUCCESS);
}
