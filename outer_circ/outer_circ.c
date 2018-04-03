/*! \file
    \brief Outer Map of the Circular Problem: main prog.

    Compute the integrals $\omega_+^*$ and $\omega_-^*$ related to the outer
    map of the circular problem.


    $Author: roldan $
    $Date: 2013-03-26 22:32:02 $
*/

#include <stdio.h>	// perror
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE
#include <assert.h>

#include <gsl/gsl_integration.h>	// gsl_integration_qags

// WARNING!!!! WE PROBABLY WANT TO USE PRTBP_DEL_CAR HERE!!!!
#include <prtbpdel.h>			// section_t

//#include <frtbpred.h>

#include <inner_circ.h>	// integrand_omega_in, iparams_omega_in

// name OF FUNCTION: omega_pm_old 
//
// WARNING! this is an old function, not used anymore.
//
// PURPOSE
// =======
// Given an energy level $H$, compute $\omega_{+,-}^*(H)$, defined as
//
//    \omega_+^*(I) = 
//       -lim_{N\to +\infty} 
//          \int_0^{14N\pi} f0(\gamma^*(s)) ds + NT_0(I),
//    \omega_-^*(I) = 
//       -lim_{N\to -\infty} 
//          \int_0^{14N\pi} f0(\gamma^*(s)) ds + NT_0(I),
//
// where $\gamma^{f,b}(s)$ is the homoclinic trajectory that starts at the
// primary homoclinic point z_2,z_1.
// 
// This is computed using numerical integration.
//
// PARAMETERS
// ==========
// mu
//    mass parameter for the RTBP
// x
//    x=(l,L,g=0,G), primary homoclinic point, on the section g=0.
//    This is actually the homoclinic point z_1 or z_2, depending if we are
//    using the outer or inner splitting.
// N
//    number of forward iterates of the poincare map (length of integration).
//    N is positive when computing $\omega_+^*$, or negative when computing
//    $\omega_-^*$.
// T0
//    shift of inner map
// omega
//    On return of this function, omega contains the value of the function.
// 
// RETURN VALUE
// ============
// Returns a non-zero error code to indicate an error and 0 to indicate
// success.
//
// NOTES
// =====
// 
// CALLS TO: integrand_omega_in

int omega_pm_old(double mu, double x[DIM], int N, double T0, double *omega)
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

   // Integrate integrand function from 0 to 14N\pi. 
   // Notice that N may be positive or negative.
   // We request a absolute error of 0 and a relative error $10^{-13}$.
   gsl_integration_qags (&F, 0, 14*N*M_PI, 0, 1.e-13, 1000, w, &result, &error);
   fprintf (stderr, "estimated error = % .3le\n", error);

   gsl_integration_workspace_free (w);

   // \omega_{+,-}^*
   *omega = -(result + N*T0);
   return 0;
}

/**
  Given an energy level \f$H\f$, compute \f$\omega_+^*(H)\f$.

  Given an energy level \f$H\f$, compute \f$\omega_+^*(H)\f$, defined as
 
  \f[ \omega_+^*(I) = 
        lim_{N\to \infty} 
           \int_0^{6N\pi} f0(\gamma^*(s)) ds + NT_0(I),
  \f]
 
  where \f$\gamma^*(s)\f$ is the homoclinic trajectory of the reduced flow that
  starts at the primary homoclinic point z.
 
  Equivalently, we compute it as
  
  \f[\omega_+^*(I) = 
        lim_{N\to +\infty} 
           \int_{-6N\pi}^0 f0(\gamma^*(s)) ds + NT_0(I),
  \f]
 
  where \f$\gamma^*(s)\f$ is the homoclinic trajectory that starts at the point 
  \f$ z^s = \mathcal{P}^{N}(z) \f$. Notice that z^s is in the STABLE manifold of
  the reduced flow.
 
  Equivalently, we split the integral from -6N\pi to 0 into N parts of size
  6\pi:
 
  \f[\omega_+^*(I) = 
        lim_{N\to +\infty} 
           \sum_{i=N,1} 
              (\int_{-6\pi}^{0} f0(\gamma_i(s)) ds + T_0(I)),
  \f]
 
  where \f$\gamma_i(s)\f$ is the homoclinic trajectory that starts at the point
  P^{-(N-i)}(z^s) = P^{i}(z). 
 
  This is computed using numerical integration.
 
  \param[in] mu 	mass parameter for the RTBP
  \param[in] x[DIM]	x=(l=0,L,g,G), point z^s, on the section l=0.

  \param[in] N
     number of iterates of the poincare map (length of integration).
     N must be positive since we are computing \f$\omega_+^*\f$.

  \param[in] T0	shift of inner map

  \param[out] omega
     On return of this function, omega contains the value of the function
     \f$\omega_+\f$.
  
  \returns
  a non-zero error code to indicate an error and 0 to indicate success.
 
  \remark We assume that point $z_s$ is on the SEC1 \f$ \{\ell=0\} \f$
  section.

 */

// NOTE: Instead of P^{-(N-i)}(z^s), we could have used
// frtbp_red(-2(N-i)\pi, z^s). They should give the same point.

int omega_pos(double mu, double x[DIM], int N, double T0, double *omega)
{
   double result, error;

   // auxiliary variables
   int i,j;
   double xi[DIM];		/* point \xi=P^{-(N-i)}(z^s) */
   double t;

   assert(N>0);

   gsl_integration_workspace * w
               = gsl_integration_workspace_alloc (1000);

   struct iparams_omega_in params;

   gsl_function F;
   F.function = &integrand_omega_in;
   F.params = &params;

   params.mu = mu;
   
   *omega = 0.0;
   for(i=N; i>=1; i--)
   {
      // $\gamma_i(s)$ is the homoclinic trajectory that starts at the
      // point \xi = P^{-(N-i)}(z^s) = P^{i}(z). 
      for(j=0;j<DIM;j++) xi[j]=x[j];

      // WARNING!!!! WE PROBABLY WANT TO USE PRTBP_DEL_CAR HERE!!!!
      if(prtbp_del_inv(mu,SEC1,(N-i)*3,xi,&t))
      {
	 fprintf(stderr, "omega_pos: error computing point P^{i}(z)\n");
	 return(1);
      }

      params.x[0] = xi[0];
      params.x[1] = xi[1];
      params.x[2] = xi[2];
      params.x[3] = xi[3];

      // Integrate integrand function from -6\pi to 0. 
      // We request a absolute error of 0 and a relative error $10^{-13}$.

      // NOTE: we can't reach accuracy of 10^{-13} when computing
      // omega_pos_f, so we lower it to 10^{-9}
      gsl_integration_qags (&F, -6*M_PI, 0, 0, 1.e-9, 1000, w, &result,
	    &error); 
      fprintf (stderr, "estimated error = % .3le\n", error);

      // \omega_+
      (*omega) = (*omega) +(result + T0);
   }

   gsl_integration_workspace_free (w);
   return 0;
}

/**
  Given an energy level \f$H\f$, compute \f$\omega_-^*(H)\f$.

  Given an energy level \f$H\f$, compute \f$\omega_-^*(H)\f$, defined as
 
  \f[ \omega_-^*(I) = 
        lim_{N\to -\infty} 
           \left( \int_0^{6N\pi} f0(\gamma^*(s)) ds + NT_0(I) \right),
  \f]
 
  where \f$\gamma^*(s)\f$ is the homoclinic trajectory of the reduced flow that
  starts at the primary homoclinic point z.
 
  Equivalently, we compute it as
  
  \f[\omega_-^*(I) = 
        lim_{N\to \infty} 
           \left( \int_{6N\pi}^0 f0(\gamma^*(s)) ds - NT_0(I) \right),
  \f]
 
  where \f$\gamma^*(s)\f$ is the homoclinic trajectory that starts at the point 
  \f$ z^u = \mathcal{P}^{-N}(z) \f$. Notice that z^u is in the UNSTABLE
  manifold of the reduced flow.
 
  Equivalently, we split the integral from 6N\pi to 0 into N parts of size
  6\pi:
 
  \f[\omega_-^*(I) = 
        lim_{N\to \infty} 
           \left( \sum_{i=N,1} 
              (\int_{6\pi}^{0} f0(\gamma_i(s)) ds - T_0(I)) \right),
  \f]
 
  where \f$\gamma_i(s)\f$ is the homoclinic trajectory that starts at the point
  P^{N-i}(z^u) = P^{-i}(z). 
 
  This is computed using numerical integration.
 
  \param[in] mu 	mass parameter for the RTBP
  \param[in] x[DIM]	x=(l=0,L,g,G), point z^u, on the section l=0.

  \param[in] N
     number of iterates of the poincare map (length of integration).
     N must be a POSITIVE integer.

  \param[in] T0	shift of inner map

  \param[out] omega
     On return of this function, omega contains the value of the function
     \f$\omega_-\f$.
  
  \returns
  a non-zero error code to indicate an error and 0 to indicate success.
 
  \remark We assume that point $z_u$ is on the SEC1 \f$ \{\ell=0\} \f$
  section.

 */

// NOTE: Instead of P^{N-i}(z^u), we could have used
// frtbp_red(2(N-i)\pi, z^u). They should give the same point.

int omega_neg(double mu, double x[DIM], int N, double T0, double *omega)
{
   double result, error;

   // auxiliary variables
   int i,j;
   double xi[DIM];		/* point \xi=P^{N-i}(z^u) */
   double t;

   assert(N>0);

   gsl_integration_workspace * w
               = gsl_integration_workspace_alloc (1000);

   struct iparams_omega_in params;

   gsl_function F;
   F.function = &integrand_omega_in;
   F.params = &params;

   params.mu = mu;
   
   *omega = 0.0;
   for(i=N; i>=1; i--)
   {
      // $\gamma_i(s)$ is the homoclinic trajectory that starts at the
      // point \xi = P^{N-i}(z^u) = P^{-i}(z). 
      for(j=0;j<DIM;j++) xi[j]=x[j];

      // WARNING!!!! WE PROBABLY WANT TO USE PRTBP_DEL_CAR HERE!!!!
      if(prtbp_del(mu,SEC1,(N-i)*3,xi,&t))
      {
	 fprintf(stderr, "omega_neg: error computing point P^{%d}(z^u)\n",
	       (N-i));
	 return(1);
      }

      params.x[0] = xi[0];
      params.x[1] = xi[1];
      params.x[2] = xi[2];
      params.x[3] = xi[3];

      // Integrate integrand function from 6\pi to 0. 
      // We request a absolute error of 0 and a relative error $10^{-13}$.

      // NOTE: we can't reach accuracy of 10^{-13} when computing
      // omega_pos_f, so we lower it to 10^{-9}
      gsl_integration_qags (&F, 6*M_PI, 0, 0, 1.e-9, 1000, w, &result,
	    &error); 
      fprintf (stderr, "estimated error = % .3le\n", error);

      // \omega_-
      (*omega) = (*omega) +(result - T0);
   }

   gsl_integration_workspace_free (w);
   return 0;
}

/**
  Outer Map of the Circular Problem: main prog.

  This program computes the integrals \f$ \omega_{+,-}^* \f$ related to the
  outer map of the circular problem.
  (See notes "Kirkwood Gaps" by Marcel).

  It reads the following input from stdin:
  - mu
     mass parameter for the RTBP

  And a sequence of lines:
  - T
     period of periodic orbit.
  - zu
     preimage of homoclinic point z (in the original flow): P^{M}(z_u) = z.
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
   double zu[DIM];	/* preimage of primary homoclinic point */
   double w_pos;	/* value of integral \omega_+^* */
   double w_neg;	/* value of integral \omega_-^* */
   double w_out;	/* value of integral \omega_{out}^* */

   double T, T0;		/* period,shift of inner map */

   /* Number of iterates to hit z from z_u: P^M(z_u)=z. 
      This will also be used as N, the number of iterates of poincare map
      along homoclinic orbit (length of integration) */
   int M;	

   // aux vars
   int i;
   double t;

   // Input mass parameter
   if(scanf("%le", &mu)<1)
   {
      perror("main: error reading input");
      exit(EXIT_FAILURE);
   }

   // Input period T, zu, number of poincare iterates M, from stdin.
   while(scanf("%le %le %le %le %le %d", &T, zu, zu+1, zu+2, zu+3, &M) == 6)
   {
      // TESTING...
      //T0 = (T-2*M_PI)/mu;
      T0 = T-2*M_PI;

      //omega_pos(mu, zu, M, T0, &w_pos);

      // Compute $\omega_-^*$, integrating along $z(s) = \gamma^*(s)$.
      // Note: since the homoclinic point is at the symmetry axis, we have
      // \omega_-^* = -\omega_+^*.
      //w_neg = -w_pos;

      // TESTING
      omega_neg(mu, zu, M-1, T0, &w_neg);
      printf("%.15e\n", w_neg);
      omega_neg(mu, zu, M, T0, &w_neg);

      //w_out = w_pos-w_neg;

      // Output result to stdout.
      printf("%.15e\n", w_neg);
      //printf("%.15e\n", w_pos);
      fflush(NULL);
   }
   exit(EXIT_SUCCESS);
}
