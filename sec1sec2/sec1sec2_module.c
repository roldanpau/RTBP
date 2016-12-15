/*! \file
  \brief Take a point $p$ on section S1 to section S2 

   Consider the RTBP. 
   Let S1 be the Poincare section {y=0, vy>0}.
   Let S2 be the Poincare section {y=0, vy<0}.
   Fix the Hamiltonian to a given value "H". 
   Using this energy condition, we can work with only two variables, $(x,p_x)$.
   The third variable $y$ is 0 since we look at the Poincare section, and the
   fourth variable $p_y$ can be obtained from the energy condition.
   These functions take a point $p$ on section S1 to section S2 by either forward
   or backward integration.

  $Author: roldan $
  $Date: 2012-12-20 11:02:51 $
  */

#include <stdio.h>	// fprintf
#include <stdlib.h>	// EXIT_FAILURE
#include <math.h>	// sqrt, fabs
#include <rtbp.h>	// DIM, rtbp
#include <hinv.h>
#include <prtbp.h>

int sec1sec2(double mu, double H, double p[2], double *ti)
{
   double x[DIM];
   x[0]=p[0];	// x
   x[1]=0.0;	// y
   x[2]=p[1]; 	// p_x

   // Compute x[3]=p_y by inverting the Hamiltonian.
   // Recall that p is on section SEC1.
   if(hinv(mu,SEC1,H,x))
   {
      fprintf(stderr, "sec1sec2: error inverting the Hamiltonian\n");
      return(1);
   }

   // Compute first intersection with section S2 (and integration time).
   // On exit, "prtbp" guarantees that x is exactly on Poincare section S2.
   if(prtbp(mu,SEC2,1,x,ti))
   {
      fprintf(stderr, "sec1sec2: error computing poincare map\n");
      return(1);
   }
   // Set the image point on section S2.
   p[0]=x[0]; p[1]=x[2];
   return(0);
}

int sec1sec2_inv(double mu, double H, double p[2], double *ti)
{
   double x[DIM];
   x[0]=p[0];	// x
   x[1]=0.0;	// y
   x[2]=p[1]; 	// p_x

   // Compute x[3]=p_y by inverting the Hamiltonian.
   // Recall that p is on section SEC1.
   if(hinv(mu,SEC1,H,x))
   {
      fprintf(stderr, "sec1sec2_inv: error inverting the Hamiltonian\n");
      return(1);
   }

   // Compute first intersection with section S2 (and integration time).
   // On exit, "prtbp_inv" guarantees that x is exactly on Poincare section S2.
   if(prtbp_inv(mu,SEC2,1,x,ti))
   {
      fprintf(stderr, "sec1sec2_inv: error computing inverse poincare map\n");
      return(1);
   }
   // Set the image point on section S2.
   p[0]=x[0]; p[1]=x[2];
   return(0);
}
