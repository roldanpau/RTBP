// ================================================
// Cartesian to Delaunay Change of Coordinates (2D)
// ================================================
// FILE:          $RCSfile: cardel2s_2d_main.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2011-07-01 08:54:47 $
//
// PURPOSE
// =======
// Obtain rotating Delaunay coordinates (l,L,g,G) from rotating cartesian
// (x,px).
//
// OBS! We DO NOT assume that all points belong to the same energy level H!
//
// OVERALL METHOD
// ==============
//
// 1. Input parameters from stdin:
// 
//    - mass parameter "mu"
//
// 2. For each energy value in input line do:
//    2.1 Input data:
//       - energy value "H"
//       - point in cartesian coordinates (x,px)
//    2.2 Obtain Delaunay coordinates
//    2.3 Output one line to stdout:
//       - energy value "H"
//       - point in Delaunay coordinates (l,L,g,G)
//
// NOTES
// =====
// Input point is on the Poincare section $\Sigma_-$. We can recover the 4
// euclidean coordinates (x,y,px,py) using the energy relation.

// Headers
#include <stdio.h>
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE
#include <rtbp.h>	// DIM
#include "cardel2_2d.h"

int main( )
{
   double mu, H; 
   double x, px;

   // auxiliary variables
   double X[2], Y[DIM];
   int status;

   // Input mass parameter
   if(scanf("%le", &mu) < 1)
   {
      perror("main: error reading input");
      exit(EXIT_FAILURE);
   }

   // Process list of points in cartesian coordinates.
   while(scanf("%le %le %le", &H, &x, &px) == 3)
   {
      // 2. Obtain Delaunay coordinates
      X[0]=x;
      X[1]=px;
      status = cardel2_2d(mu,H,X,Y);
      if(status)
      {
	 fprintf(stderr, "main: error obtaining Delaunay coordinates\n");
	 exit(EXIT_FAILURE);
      }

      // 3. Output point in Delaunay to stdout.
      if(printf("%e %.15e %.15e %.15e %.15e\n", H, Y[0], Y[1], Y[2], Y[3])<0)
      {
	 perror("main: error writting output");
	 exit(EXIT_FAILURE);
      }
   }
   exit(EXIT_SUCCESS);
}
