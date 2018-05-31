/*! \file
    \brief Outer Map of the Circular Problem for Stochastic paper: test prog.

    Test convergence of the integral defining $\omega_-$ as a function of the
    number of poincare iterates $N$ in the following way:
    For each $N$, compute the difference $|\int_0^N ... - \int_0^^{N-1}|$.

    This difference should go to zero as $N\to +\infty$, and ideally it should
    be much smaller than the size of $\omega_-$ itself.

    $Author: roldan $
    $Date: 2013-03-26 22:32:02 $
*/

#include <stdio.h>	// perror
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE
#include <string.h> // strcmp
#include <assert.h>
#include <math.h>   // M_PI

#include <section.h>
#include "outer_circ_stoch_module.h"

/**
  Outer Map of the Circular Problem: test prog.

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
