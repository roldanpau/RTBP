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
