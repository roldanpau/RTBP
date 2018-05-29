/*! \file rtbpred.c
    \brief Reduced Restricted Three Body Problem equations.

    $Author: roldan $
    $Date: 2013-03-26 22:18:24 $
*/

#include <math.h>	        // fabs
#include <gsl/gsl_errno.h>	// GSL_SUCCESS
#include <rtbpdel.h>		// ERR_COLLISION, rtbp_del

/**
  Computes the vectorfield of the reduced RTBP problem.

  The function computes the vectorfield of the reduced RTBP at a point
  named "x".
  Specifically, we use the equations of motion of the reduced planar
  circular RTBP in Delaunay coordinates (l,L,g,G,t,I), in the rotating
  (synodic) coordinate system.

  \param[in] s
  adimensional time at which the vectorfield is evaluated. Since this is
  an autonomous ODE (does not depend on time), this parameter is not used.

  \param[in] x[6]
  point in phase space, 6 coordinates: (l,L,g,G,t,I).

  \param[out] y[6]
  vectorfield at (s,x), 6 coordinates: d/ds(l,L,g,G,t,I).

  \param[in] params
  pointer to the parameter of the system: the mass ratio "mu".

  \returns
  status code of the function (success/error).

  \retval GSL_SUCCESS success.
  \retval 1 Error computing non-reduced vector field. 

  \author Marcel Guardia
  \author Pau Roldan

  \remark
  See Marcel's notes "Inner and outer dynamics" and "Kirkwood Gaps".
 */

int rtbp_red_l(double s, const double *x, double *y, void *params)
{
   // auxiliary variables
   double dl;			// dl/dt

   // In the call to rtbp_del, only the 4 first components of x,y are used.
   // WARNING: Should we call rtbp_del with $s$ as independent variable, or
   // with $t$?? This does not matter for the circular case, since it is
   // autonomous, but probably it matters for the elliptic case.
   if(rtbp_del(s,x,y,params))
   {
      fprintf(stderr, "rtbp_red: error computing non-reduced vector field\n");
      return 1;
   }
   // Reduced vector field is simply the non-reduced one scaled by the $\dot
   // l$ component.
   dl = y[0];
   if(fabs(dl)<1.e-15)
   {
       fprintf(stderr, "rtbp_red: Singularity of the vectorfield!\n");
       return 2;
   }

   y[0] = 1; 		// dl/ds
   y[1] /= dl;		// dL/ds
   y[2] /= dl; 		// dg/ds
   y[3] /= dl;		// dG/ds
   y[4] = 1.0/dl;	// dt/ds
   y[5] = 0;		// dI/ds

   return GSL_SUCCESS;
}

int rtbp_red_g(double s, const double *x, double *y, void *params)
{
   // auxiliary variables
   double dg;			// dg/dt

   // In the call to rtbp_del, only the 4 first components of x,y are used.
   // WARNING: Should we call rtbp_del with $s$ as independent variable, or
   // with $t$?? This does not matter for the circular case, since it is
   // autonomous, but probably it matters for the elliptic case.
   if(rtbp_del(s,x,y,params))
   {
      fprintf(stderr, "rtbp_red: error computing non-reduced vector field\n");
      return 1;
   }
   // Reduced vector field is simply the non-reduced one scaled by the $\dot
   // g$ component.
   dg = y[2];
   if(fabs(dg)<1.e-15)
   {
       fprintf(stderr, "rtbp_red: Singularity of the vectorfield!\n");
       return 2;
   }

   // DEBUG:
   fprintf(stderr, "dl=%e, dg=%e\n", y[0], y[2]);

   y[0] /= dg; 		// dl/ds
   y[1] /= dg;		// dL/ds
   y[2] = 1; 		// dg/ds
   y[3] /= dg;		// dG/ds
   y[4] = 1.0/dg;	// dt/ds
   y[5] = 0;		// dI/ds

   return GSL_SUCCESS;
}
