/*! \file
    \brief Flow of the Restricted Three Body Problem

    $Author: roldan $
    $Date: 2013-03-26 22:17:52 $
*/

#include <taylor2d.h>	// taylor_step_rtbp2d
#include <taylor2dv.h>	// taylor_step_rtbp2dv
#include <rtbp.h>	// DIM
#include "frtbp.h"	// DIMV

double mu;	///< global variable, needed in taylor_step_rtbp2d and taylor_step_rtbp2dv

// NOTES
// =====
// Notice that the taylor integrator was written using the opposite
// convention. To avoid re-writting the taylor integrator, we simply call
// taylor with parameter mu=1-mu.

int dfrtbp(double mu_loc, double t1, double x[DIM], double dphi[DIMV])
{
   /* (log10) absolute error for local error control */
   double log10_eps_abs = -16;	

   /* (log10) relative error for local error control */
   double log10_eps_rel = -16;	

   double t = 0.0;
   double h;		/* step size */
   double xvar[DIM+DIMV];	/* point + derivative */
   int status, i, nt;

   // Set global variable "mu". This needs to defined before calling
   // taylor_step_rtbp2dv
   mu=1.0-mu_loc;

   // Initial condition
   for(i=0; i<DIM; i++)
      xvar[i] = x[i];

   // Initialize derivative to the identity.
  xvar[ 4]=1; xvar[ 5]=0; xvar[ 6]=0; xvar[ 7]=0;
  xvar[ 8]=0; xvar[ 9]=1; xvar[10]=0; xvar[11]=0;
  xvar[12]=0; xvar[13]=0; xvar[14]=1; xvar[15]=0;
  xvar[16]=0; xvar[17]=0; xvar[18]=0; xvar[19]=1;

   // Integrate trajectory numerically.
   // The time $t1$ may be positive or negative, allowing for forward or
   // backward integration.
   if(t1>=0)
   {
      // Forward integration
      while (t<t1)
      {
	    status =
        taylor_step_rtbp2dv(&t,xvar,1,1,log10_eps_abs,log10_eps_rel,&t1,&h,&nt,NULL);
      }
   }
   else // t1<0
   {
      // Backwards integration
      while (t>t1)
	    status =
        taylor_step_rtbp2dv(&t,xvar,-1,1,log10_eps_abs,log10_eps_rel,&t1,&h,&nt,NULL);
   }

   // Save derivative to matrix "dphi"
   for(i=0; i<DIMV; i++)
      dphi[i] = xvar[DIM+i];
   return(0);
}

// NOTES
// =====
// Notice that the taylor integrator was written using the opposite
// convention. To avoid re-writting the taylor integrator, we simply call
// taylor with parameter mu=1-mu.

int frtbp(double mu_loc, double t1, double x[DIM])
{
   double log10_eps_abs = -16;	/* (log10) absolute error for local error control */
   double log10_eps_rel = -16;	/* (log10) relative error for local error control */

   double t = 0.0;
   double h;		/* step size */
   int status, nt;

   // Set global variable "mu". This needs to defined before calling
   // taylor_step_rtbp2d
   mu=1.0-mu_loc;

   // Integrate trajectory numerically.
   // The time $t1$ may be positive or negative, allowing for forward or
   // backward integration.
   if(t1>=0)
   {
      // Forward integration
      while (t<t1)
      {
	    status =
        taylor_step_rtbp2d(&t,x,1,1,log10_eps_abs,log10_eps_rel,&t1,&h,&nt,NULL);
	 // What if there is a singularity in the vectorfield?? This is not
	 // contemplated in the errorcode status??
      }
   }
   else // t1<0
   {
      // Backwards integration
      while (t>t1)
	    status =
        taylor_step_rtbp2d(&t,x,-1,1,log10_eps_abs,log10_eps_rel,&t1,&h,&nt,NULL);
   }
   return(0);
}
