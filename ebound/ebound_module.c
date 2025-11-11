/*! \file ebound_module.c
    \brief Bounds for \f$ \mathcal{E}_1^1(X,t) \f$ and \f$ E_2(X) \f$.
    
    Compute bound of 
    \f$ \mathcal{E}_1^1(X,t) = E(J,L_j(t,J)) - E(J,L_0) \f$ for
    all $\fJ\in[J_-, J_+]\f$, where the function E is defined as
    \f[ E(J,L) = \sqrt{1-\frac{(J+\frac{1}{2L^2})^2}{L^2}}. \f]

    Also compute bound of 
    \f$ E_2(X) =  \sqrt{1-\frac{G^2}{L^2}} - E(J,L) \f$ for
    all $\fJ\in[J_-, J_+]\f$, where the function E is defined above.

    \note In principle, \f$ E(J,L_j(t,J)) \f$ varies along the infinite-time
    homoclinic. However the homoclinic tends to the periodic orbit, so once we
    are close to the periodic orbit, it is okay to disregard both ``tails'' of
    the homoclinic.
*/

#include <stdio.h>	// perror
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE

#include <utils_module.h>           // dblcpy
#include <frtbp.h>
#include <cardel.h>
#include <math.h>           // sqrt, pow, cbrt

const double SHORT_TIME=0.01;       ///< integration "step" for frtbp

double E(double J, double L)
{
   double aux = J+1.0/(2*L*L);
   return sqrt(1.0 - aux*aux / (L*L));
}

int ebound(double mu, double J, double x[DIM], double t, double *E11_bound,
double *E2_bound) 
{
    double L0 = 1.0/cbrt(3);    ///< resonant value \f$ L_0 \f$

   // auxiliary variables
   int i, status;
   double xt[DIM];		    /* point x(t) in Cartesian */
   double xt_del[DIM];		/* point x(t) in Delaunay*/
   double dt, L, G, E11, E2;

   dt = (t>0 ? SHORT_TIME : -SHORT_TIME);

   /* Work with local copy to avoid modifying original x */
   dblcpy(xt, x, DIM);
   
   *E11_bound = 0.0;
   *E2_bound = 0.0;
   for(i=0; i<(t/dt); i++)
   {
	   status = frtbp(mu,dt,xt);
	   if(status)
	   {
         fprintf(stderr, "ebound: integration error during iteration %d\n",
               i);
         return(1);
	   }
	   cardel(xt,xt_del);
       L = xt_del[1];
       G = xt_del[3];

       /* Update bounds E11 and E2*/

       E11 = E(J,L)-E(J,L0);
       if(fabs(E11) > *E11_bound) *E11_bound=fabs(E11);

       E2 = sqrt(1.0 - G*G / (L*L)) - E(J,L);
       if(fabs(E2) > *E2_bound) *E2_bound=fabs(E2);
   }
   return 0;
}
