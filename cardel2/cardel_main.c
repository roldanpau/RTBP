// ===========================================
// Cartesian to Delaunay Change of Coordinates
// ===========================================
// FILE:          $RCSfile: cardel_main.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2011-03-31 13:07:57 $
//
// PURPOSE
// =======
// Obtain rotating Delaunay coordinates (l,L,g,G) from rotating cartesian
// (x,y,px,py).
//
// OVERALL METHOD
// ==============
//
// 1. Input parameters from stdin:
// 
//    - mass parameter "mu"
//    - energy value "H"
//    - point in cartesian coordinates (x,y,px,py)
//
// 2. Obtain Delaunay coordinates
// 3. Output point in Delaunay to stdout.
//
// NOTES
// =====

// Headers
#include <stdio.h>
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE
#include <rtbp.h>	// DIM
#include <hinv.h>
#include "cardel.h"

int main( )
{
   double mu;		// actually, mu is not used
   double H; 
   double x, y, px, py;

   // auxiliary variables
   double X[DIM], Y[DIM];
   int status;

   // Input mass parameter, energy value, momentum estimate.
   if(scanf("%le %le", &mu, &H) < 2)
   {
      perror("main: error reading input");
      exit(EXIT_FAILURE);
   }

   // Process list of points in cartesian coordinates.
   while(scanf("%le %le %le %le", &x, &y, &px, &py) == 4)
   {
      X[0] = x;
      X[1] = y;
      X[2] = px;
      X[3] = py;

      // 2. Obtain Delaunay coordinates
      cardel(X,Y);

      // 3. Output point in Delaunay to stdout.
      if(printf("%.15e %.15e %.15e %.15e\n", Y[0], Y[1], Y[2], Y[3])<0)
      {
	 perror("main: error writting output");
	 exit(EXIT_FAILURE);
      }
   }
   exit(EXIT_SUCCESS);
}
