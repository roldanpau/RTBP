/*! \file Lbound.c
    \brief Bound for \f$L_{hom}(t,J)\f$: main program
    
    Compute bound \f$|L_{hom}(t,J) - L_0| \leq C\mu\f$ for all $\fJ\in[J_-,
    J_+]\f$.

    \note In principle, L(t) varies along the infinite-time homoclinic. However
    the homoclinic tends to the periodic orbit, so once we are close to the
    periodic orbit, it is okay to disregard both ``tails'' of the homoclinic.
*/

#include <stdio.h>	// perror
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE
#include <string.h> // strcmp
#include <assert.h>
#include <math.h>   // M_PI

#include <frtbp.h>
#include <approxint.h>	// stability_t
#include "Lbound_module.h"

/** 
Bound for \f$L_{hom}(t,J)\f$: main program

Compute an upper bound for \f$L_{hom}(t)\f$ for all $\fJ\in[J_-, J_+]\f$.

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
     (Cartesian coordinates). Note: WE USE THIS SAME VARIABLE FOR z_s.
  - t 
     integration time to reach homoclinic point z from z_u

  For each input line, it outputs result to stdout:
  - bound
 
 */
 
int main( )
{
   double mu;

   // "stability" flag specifies wheather we want to use the unstable branch
   // (=0) or stable branch (=1) of the manifold
   int stability;

   double zu[DIM];	/* preimage of primary homoclinic point */
   double bound;	/* bound \f$ |L_{hom}(t)| \f$ */
   double T;	    /* period */
   double t;	    /* integration time from z_u to z */

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
	   if(st==UNSTABLE)
	   {
		   // Instead of fetching zs from intersecs_st_SECg_br1.res, 
		   // we set it to the symmetric point of zu. 
		   // This is done so that omega_neg and omega_pos are actually symmetric.
		   /*
		   dblcpy(zs_car,zu_car,DIM);
		   zs_car[1] = -zu_car[1];
		   zs_car[2] = -zu_car[2];
		   */

		  // Compute Lbound, integrating along $z(s) = \gamma^*(s)$.
		  Lbound(mu, zu, t, &bound);

		  // Output result to stdout.
		  printf("%.15e\n", bound);
		  fflush(NULL);
	   }
	   else // st==STABLE
	   {
           // t should be negative
           assert(t<0);

		   // Instead of fetching zs from intersecs_st_SECg_br1.res, 
		   // we set it to the symmetric point of zu. 
		   // This is done so that omega_neg and omega_pos are actually symmetric.
		   /*
		   dblcpy(zs_car,zu_car,DIM);
		   zs_car[1] = -zu_car[1];
		   zs_car[2] = -zu_car[2];
		   */

		  // Compute Lbound, integrating along $z(s) = \gamma^*(s)$.
		  Lbound(mu, zu, t, &bound);   

		  // Output result to stdout.
		  printf("%.15e\n", bound);
		  fflush(NULL);
	   }
   }
   exit(EXIT_SUCCESS);
}
