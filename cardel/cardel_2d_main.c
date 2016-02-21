/*! \file
    \brief Cartesian to Delaunay Change of Coordinates (2D): main prog.

    Obtain rotating Delaunay coordinates (l,L,g,G) from rotating cartesian
    (x,px).

    $Author: roldan $
    $Date: 2013-03-26 22:15:14 $
*/

// Headers
#include <stdio.h>
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE
#include <rtbp.h>	// DIM
#include "cardel_2d.h"

/**
  Cartesian to Delaunay Change of Coordinates (2D): main prog.

  Obtain rotating Delaunay coordinates (l,L,g,G) from rotating cartesian
  (x,px).
  Input point is on the Poincare section (y=0). We can recover the 4
  euclidean coordinates (x,y,px,py) using the energy relation.

  \remark
  We assume that all points belong to the same energy level H!

  1. Input parameters from stdin:
    - mass parameter "mu"
    - energy value "H"

  2. For each point in cartesian coordinates (x,px) read from stdin, do
    2.1 Obtain Delaunay coordinates
    2.2 Output point in Delaunay to stdout.

 */
int main( )
{
   double mu, H; 
   double x, px;

   // auxiliary variables
   double X[2], Y[DIM];
   int status;

   // Input mass parameter, energy value.
   if(scanf("%le %le", &mu, &H) < 2)
   {
      perror("main: error reading input");
      exit(EXIT_FAILURE);
   }

   // Process list of points in cartesian coordinates.
   while(scanf("%le %le", &x, &px) == 2)
   {
      // 2. Obtain Delaunay coordinates
      X[0]=x;
      X[1]=px;
      status = cardel_2d(mu,H,X,Y);
      if(status)
      {
	 fprintf(stderr, "main: error obtaining Delaunay coordinates\n");
	 exit(EXIT_FAILURE);
      }

      // 3. Output point in Delaunay to stdout.
      if(printf("%.15e %.15e %.15e %.15e\n", Y[0], Y[1], Y[2], Y[3])<0)
      {
	 perror("main: error writting output");
	 exit(EXIT_FAILURE);
      }
   }
   exit(EXIT_SUCCESS);
}
