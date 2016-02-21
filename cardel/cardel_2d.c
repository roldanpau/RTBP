// Headers
#include <stdio.h>	// fprintf
#include <rtbp.h>	// DIM
#include <hinv.h>
#include "cardel.h"

int cardel_2d(double mu, double H, double z[2], double Y[DIM])
{
   // Point is on Poincare section {y=0}
   double y=0;

   // auxiliary variables
   double X[DIM];
   int status;

   // Recover py using energy condition
   X[0] = z[0];
   X[1] = y;
   X[2] = z[1];
   status=hinv(mu,SEC2,H,X);
   if(status)
   {
      fprintf(stderr, \
	    "cardel_2d: error inverting Hamiltonian equation\n");
      return(1);
   }

   // Obtain Delaunay coordinates
   cardel(X,Y);
   return(0);
}
