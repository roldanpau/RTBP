/*! \file
    \brief Cartesian to Delaunay Change of Coordinates: main prog.

    Obtain rotating Delaunay coordinates (l,L,g,G) from rotating cartesian
    (x,y,px,py).

    $Author: roldan $
    $Date: 2013-03-26 22:15:14 $
*/

// Headers
#include <stdio.h>
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE
#include <rtbp.h>	// DIM
#include "cardel.h"

int main( )
{
   double x, y, px, py;

   // auxiliary variables
   double X[DIM], Y[DIM];
   int status;

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
