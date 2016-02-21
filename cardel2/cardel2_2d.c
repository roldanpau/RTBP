// ================================================
// Cartesian to Delaunay Change of Coordinates (2D)
// ================================================
// FILE:          $RCSfile: cardel2_2d.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2011-03-31 13:07:57 $
//
// PURPOSE:
//
// NOTES:
//
// OVERALL METHOD:
//
// FUNCTIONS
// =========
//
// cardel2_2d
// ---------
// Obtain rotating Delaunay coordinates (l,L,g,G) from rotating cartesian
// (x,px).

// Headers
#include <stdio.h>	// fprintf
#include <rtbp.h>	// DIM
#include <hinv2.h>
#include "cardel.h"

// name OF FUNCTION: cardel2_2d
//
// PURPOSE
// =======
// Obtain rotating Delaunay coordinates (l,L,g,G) from rotating cartesian
// (x,px).
//
// PARAMETERS
// ==========
// mu
//    mass parameter for the RTBP
// H
//    energy value
// z=(x,px)
//    rotating cartesian coordinates (input)
// Y=(l,L,g,G)
//    rotating Delaunay coordinates (output)
// 
// RETURN VALUE
// ============
// Returns a non-zero error code to indicate an error and 0 to indicate
// success.
//
// NOTES
// =====
// Input point is on the Poincare section (y=0). We can recover the 4
// euclidean coordinates (x,y,px,py) using the energy relation.
// 
// CALLS TO: hinv2, cardel
int cardel2_2d(double mu, double H, double z[2], double Y[DIM])
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
   status=hinv2(mu,H,X);
   if(status)
   {
      fprintf(stderr, \
	    "cardel2_2d: error inverting Hamiltonian equation\n");
      return(1);
   }

   // Obtain Delaunay coordinates
   cardel(X,Y);
   return(0);
}
