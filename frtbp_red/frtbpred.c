/*! \file
    \brief Flow of the Reduced Restricted Three Body Problem

    $Author: roldan $
    $Date: 2013-03-26 22:18:24 $
*/

#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>
#include "rtbpred.h"	// DIMRED, rtbp_red

double EPS_ABS=1.e-16;     /* absolute error for local error control */
double EPS_REL=0.0;        /* relative error for local error control */

int frtbp_red_l(double mu, double s1, double x[DIMRED])
{
   double t = 0.0;
   double h;    /* step size */
   int status;

   // Embedded Runge-Kutta Prince-Dormand (8,9) method.
   const gsl_odeiv_step_type *T = gsl_odeiv_step_rk8pd;
   gsl_odeiv_step *s = gsl_odeiv_step_alloc(T,DIMRED);

   // Control to determine optimal step size: keep the local error on each
   // step within an absolute error of EPS_ABS and relative error of EPS_REL
   // with respect to the solution.
   gsl_odeiv_control *c = gsl_odeiv_control_y_new(EPS_ABS,EPS_REL);
   gsl_odeiv_evolve *e = gsl_odeiv_evolve_alloc(DIMRED);

   // define system of equations (NULL = we don't provide the jacobian)
   gsl_odeiv_system sys = {rtbp_red_l,NULL,DIMRED,&mu};

   // Integrate trajectory numerically.
   // The time $s1$ may be positive or negative, allowing for forward or
   // backward integration.
   if(s1>=0)
   {
      // Forward integration
      h = 1.e-6;
      while (t<s1)
      {
         status = gsl_odeiv_evolve_apply(e,c,s,&sys,&t,s1,&h,x);
         if (status != GSL_SUCCESS)
         {
            fprintf(stderr, "frtbp_red: error integrating trajectory");
            return(1);
         }
      }
   }
   else // s1<0
   {
      // Backward integration
      h = -1.e-6;
      while (t>s1)
      {
         status = gsl_odeiv_evolve_apply(e,c,s,&sys,&t,s1,&h,x);
         if (status != GSL_SUCCESS)
         {
            fprintf(stderr, "frtbp_red: error integrating trajectory");
            return(1);
         }
      }
   }
   gsl_odeiv_evolve_free(e);
   return(0);
}

int frtbp_red_g(double mu, double s1, double x[DIMRED])
{
   double t = 0.0;
   double h;    /* step size */
   int status;

   // Embedded Runge-Kutta Prince-Dormand (8,9) method.
   const gsl_odeiv_step_type *T = gsl_odeiv_step_rk8pd;
   gsl_odeiv_step *s = gsl_odeiv_step_alloc(T,DIMRED);

   // Control to determine optimal step size: keep the local error on each
   // step within an absolute error of EPS_ABS and relative error of EPS_REL
   // with respect to the solution.
   gsl_odeiv_control *c = gsl_odeiv_control_y_new(EPS_ABS,EPS_REL);
   gsl_odeiv_evolve *e = gsl_odeiv_evolve_alloc(DIMRED);

   // define system of equations (NULL = we don't provide the jacobian)
   gsl_odeiv_system sys = {rtbp_red_g,NULL,DIMRED,&mu};

   // Integrate trajectory numerically.
   // The time $s1$ may be positive or negative, allowing for forward or
   // backward integration.
   if(s1>=0)
   {
      // Forward integration
      h = 1.e-6;
      while (t<s1)
      {
         status = gsl_odeiv_evolve_apply(e,c,s,&sys,&t,s1,&h,x);
         if (status != GSL_SUCCESS)
         {
            fprintf(stderr, "frtbp_red: error integrating trajectory");
            return(1);
         }
      }
   }
   else // s1<0
   {
      // Backward integration
      h = -1.e-6;
      while (t>s1)
      {
         status = gsl_odeiv_evolve_apply(e,c,s,&sys,&t,s1,&h,x);
         if (status != GSL_SUCCESS)
         {
            fprintf(stderr, "frtbp_red: error integrating trajectory");
            return(1);
         }
      }
   }
   gsl_odeiv_evolve_free(e);
   return(0);
}
