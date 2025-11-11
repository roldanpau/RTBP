/*! \file ebound.c
    \brief Bounds for \f$ \mathcal{E}_1^1(X,t) \f$ and \f$ E_2(X) \f$.

    \note In principle, \f$ E(J,L_j(t,J)) \f$ varies along the infinite-time
    homoclinic. However the homoclinic tends to the periodic orbit, so once we
    are close to the periodic orbit, it is okay to disregard both ``tails'' of
    the homoclinic.
*/

#include <stdio.h>	// perror
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE
#include <string.h> // strcmp
#include <assert.h>
#include <math.h>   // M_PI

#include <frtbp.h>
#include <approxint.h>	// stability_t
#include "ebound_module.h"

/** 
Bounds for \f$ \mathcal{E}_1^1(X,t) \f$ and \f$ E_2(X) \f$: main program.

Compute upper bounds for \f$ \mathcal{E}_1^1(X,t) \f$ and \f$ E_2(X) \f$ for all
$\fJ\in[J_-, J_+]\f$.

  It reads the following input from stdin:
  - mu
     mass parameter for the RTBP
  - stability
     unst/st flag

  And a sequence of lines:
  - H
     energy level: periodic (and homoclinic) orbit has energy H
  - T
     period of periodic orbit.
  - zu
     preimage of homoclinic point z (in the original flow): P^{M}(z_u) = z.
     (Cartesian coordinates). Note: WE USE THIS SAME VARIABLE FOR z_s.
  - t 
     integration time to reach homoclinic point z from z_u

  For each input line, it outputs result to stdout:
  - H
  - E11 bound
  - E2 bound
 */
 
int main( )
{
   double mu;

   // "stability" flag specifies wheather we want to use the unstable branch
   // (=0) or stable branch (=1) of the manifold
   int stability;

   double zu[DIM];	/* preimage of primary homoclinic point */

   double E11_bound;	/* bound for \f$ \mathcal{E}_1^1(X,t) \f$ */
   double E2_bound;	    /* bound for \f$ E_2(X) \f$ */

   double H;	    /* energy */
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

   // Input energy H, period T, zu, time to reach hom. pt. z, from stdin.
   while(scanf("%le %le %le %le %le %le %le", &H, &T, zu, zu+1, zu+2, zu+3, &t) == 7)
   {
	   if(st==UNSTABLE)
           assert(t>0);
       else
           assert(t<0);

       // Instead of fetching zs from intersecs_st_SECg_br1.res, 
       // we set it to the symmetric point of zu. 
       // This is done so that omega_neg and omega_pos are actually symmetric.
       /*
       dblcpy(zs_car,zu_car,DIM);
       zs_car[1] = -zu_car[1];
       zs_car[2] = -zu_car[2];
       */

       // Compute bounds, integrating along $z(s) = \gamma^*(s)$.
       ebound(mu, H, zu, t, &E11_bound, &E2_bound);

       // Output result to stdout.
       printf("%.15e %.15e %24.16e\n", H, E11_bound, E2_bound);
       fflush(NULL);
   }
   exit(EXIT_SUCCESS);
}
