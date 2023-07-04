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
#include <math.h>   // M_PI

#include <frtbp.h>
#include <approxint.h>	// stability_t
#include "outer_circ_stoch_module.h"

/**
  Outer Map of the Circular Problem: main prog.

  This program computes the integrals \f$ \omega_{+,-}^j \f$ related to the
  outer map of the circular problem.
  (See the new paper ``Stochastic Diffusion''.)

  It reads the following input from stdin:
  - mu
     mass parameter for the RTBP
  - stability
     unst/st flag

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
  - omega_pos
 
 */
 
//  NOTE: Due to symmetry, we have the following relation: 
//     \omega_- = - \omega_+, 
//  so it is enough to compute one of them.
 
int main( )
{
   double mu;

   // "stability" flag specifies wheather we want to use the unstable branch
   // (=0) or stable branch (=1) of the manifold
   int stability;

   double zu[DIM];	    /* preimage of primary homoclinic point */

   //double zs[DIM];	    /* preimage of primary homoclinic point */
   //double zs_car[DIM];	/* preimage of primary homoclinic point (Cartesian) */

   double w_pos;	/* value of integral \omega_+^* */
   double w_neg;	/* value of integral \omega_-^* */
   double w_out;	/* value of integral \omega_{out}^* */

   double T, T0;		/* period,shift of inner map */

   /* Number of iterates to hit z from z_u: P^M(z_u)=z. 
      This will also be used as N, the number of iterates of poincare map
      along homoclinic orbit (length of integration) */
   int M;	

   double t;	/* integration time from z_u to z */

   // auxiliary vars
   int status;
   stability_t st;
   double t_aux;

   // Input parameters from stdin.
   if(scanf("%le %d", &mu, &stability)<2)
   {
      perror("main: error reading input");
      exit(EXIT_FAILURE);
   }

   st = (stability==0 ? UNSTABLE : STABLE);

   // Input period T, zu, time to reach hom. pt. z, from stdin.
   while(scanf("%le %le %le %le %le %le", &T, zu, zu+1, zu+2, zu+3, &t) == 6)
   {
	  // TESTING...
	  //T0 = (T-2*M_PI)/mu;
	  T0 = T-2*M_PI;

	   if(st==UNSTABLE)
	   {
		   // Instead of fetching zs from intersecs_st_SECg_br1.res, 
		   // we set it to the symmetric point of zu. 
		   // This is done so that omega_neg and omega_pos are actually symmetric.
		   /*
		   dblcpy(zs,zu,DIM);
		   zs[0] = -zu[0];
		   zs[2] = -zu[2];

		   dblcpy(zs_car,zu_car,DIM);
		   zs_car[1] = -zu_car[1];
		   zs_car[2] = -zu_car[2];
		   */

		  // Translate zu to an exact preimage of z by the Poincare map satisfying
		  //    flow_red_g(2\pi*M, zu) = z.
		  t_aux = fmod(t,T);
		  M = t/T;

		   status=frtbp(mu,t_aux,zu);
		   if(status)
		   {
			  fprintf(stderr, "main: error integrating trajectory");
			  exit(EXIT_FAILURE);
		   }

		  // Compute $\omega_-^*$, integrating along $z(s) = \gamma^*(s)$.
		  // Note: since the homoclinic point is at the symmetry axis, we have
		  // \omega_-^* = -\omega_+^*.
		  omega_neg_stoch(mu, zu, M, T0, &w_neg);

		  // Compute $\omega_+^*$, integrating along $z(s) = \gamma^*(s)$.
		  // Note: since the homoclinic point is at the symmetry axis, we should
		  // have \omega_-^* = -\omega_+^*.

		  //omega_pos_stoch(mu, zs, M, T0, &w_pos);

		  //w_out = w_pos-w_neg;

		  // Output result to stdout.
		  printf("%.15e\n", w_neg);
		  //printf(" %.15e\n", w_pos);
		  fflush(NULL);
	   }
	   else // st==STABLE
	   {
		   // Instead of fetching zs from intersecs_st_SECg_br1.res, 
		   // we set it to the symmetric point of zu. 
		   // This is done so that omega_neg and omega_pos are actually symmetric.
		   /*
		   dblcpy(zs,zu,DIM);
		   zs[0] = -zu[0];
		   zs[2] = -zu[2];

		   dblcpy(zs_car,zu_car,DIM);
		   zs_car[1] = -zu_car[1];
		   zs_car[2] = -zu_car[2];
		   */

		  // Translate zs to an exact preimage of z by the Poincare map
		  // satisfying flow_red_g(-2\pi*M, zs) = z.
		  t_aux = fmod(t,T);	// t_aux should be negative
		  M = t/T;				// M should be negative

		   status=frtbp(mu,t_aux,zu);
		   if(status)
		   {
			  fprintf(stderr, "main: error integrating trajectory");
			  exit(EXIT_FAILURE);
		   }

		  // Compute $\omega_+^*$, integrating along $z(s) = \gamma^*(s)$.
		  // Note: since the homoclinic point is at the symmetry axis, we have
		  // \omega_-^* = -\omega_+^*.
		  omega_pos_stoch(mu, zu, -M, T0, &w_pos);

		  // Compute $\omega_+^*$, integrating along $z(s) = \gamma^*(s)$.
		  // Note: since the homoclinic point is at the symmetry axis, we should
		  // have \omega_-^* = -\omega_+^*.

		  //omega_pos_stoch(mu, zs, M, T0, &w_pos);

		  //w_out = w_pos-w_neg;

		  // Output result to stdout.
		  printf("%.15e\n", w_pos);
		  //printf(" %.15e\n", w_pos);
		  fflush(NULL);
	   }
   }
   exit(EXIT_SUCCESS);
}
