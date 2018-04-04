/*! \file
  \brief Take a point $p$ on section S1 to section S2 

   Consider the RTBP. 
   Let S1 be the Poincare section {y=0, vy>0}.
   Let S2 be the Poincare section {y=0, vy<0}.
   These functions take a point $p$ on section S1 to section S2 by either forward
   or backward integration.

  $Author: roldan $
  $Date: 2012-12-20 11:02:51 $
  */

#include <stdio.h>	// fprintf
#include <stdlib.h>	// EXIT_FAILURE
#include <math.h>	// sqrt, fabs
#include <rtbp.h>	// DIM, rtbp
#include <prtbp.h>

int sec1sec2(double mu, double p[DIM], double *ti)
{
   // Compute first intersection with section S2 (and integration time).
   // On exit, "prtbp" guarantees that x is exactly on Poincare section S2.
   if(prtbp(mu,SEC2,1,p,ti))
   {
      fprintf(stderr, "sec1sec2: error computing poincare map\n");
      return(1);
   }
   return(0);
}

int sec1sec2_inv(double mu, double p[DIM], double *ti)
{
   // Compute first intersection with section S2 (and integration time).
   // On exit, "prtbp_inv" guarantees that x is exactly on Poincare section S2.
   if(prtbp_inv(mu,SEC2,1,p,ti))
   {
      fprintf(stderr, "sec1sec2_inv: error computing inverse poincare map\n");
      return(1);
   }
   return(0);
}
