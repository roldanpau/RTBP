/*! \file Lbound_module.c
    \brief Bound for \f$L_{hom}(t,J)\f$
    
    Compute bound \f$|L_{hom}(t,J) - L_0| \leq C\mu\f$ for all $\fJ\in[J_-,
    J_+]\f$.

    \note In principle, L(t) varies along the infinite-time homoclinic. However
    the homoclinic tends to the periodic orbit, so once we are close to the
    periodic orbit, it is okay to disregard both ``tails'' of the homoclinic.
*/

#include <stdio.h>	// perror
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE

#include <utils_module.h>           // dblcpy
#include <frtbp.h>
#include <cardel.h>
#include <math.h>           // cbrt

const double SHORT_TIME=0.01;       ///< integration "step" for frtbp

int Lbound(double mu, double x[DIM], double t, double *bound) 
{
    double L0 = 1.0/cbrt(3);    ///< resonant value \f$ L_0 \f$

   // auxiliary variables
   int i, status;
   double xt[DIM];		    /* point x(t) in Cartesian */
   double xt_del[DIM];		/* point x(t) in Delaunay*/
   double dt, L;

   dt = (t>0 ? SHORT_TIME : -SHORT_TIME);

   /* Work with local copy to avoid modifying original x */
   dblcpy(xt, x, DIM);
   
   *bound = 0.0;
   for(i=0; i<(t/dt); i++)
   {
	   status = frtbp(mu,dt,xt);
	   if(status)
	   {
         fprintf(stderr, "Lbound_unst: integration error during iteration %d\n",
               i);
         return(1);
	   }
	   cardel(xt,xt_del);

       /* Update L bound */
       L = xt_del[1];
       if(fabs(L-L0) > *bound) *bound=L;
   }
   return 0;
}
