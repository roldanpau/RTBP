/*! \file
    \brief Flow of the RTBP in Delaunay coords: main prog.

    $Author: roldan $
    $Date: 2013-03-26 22:18:08 $
*/

#include <stdio.h>
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE
#include "frtbpdel.h"

/**
  Flow of the RTBP in Delaunay coords: main prog.

  Compute the flow $\phi(t,x)$ of the RTBP for a given time and initial
  condition. 
  The time $t$ may be positive or negative, allowing for forward or backward
  integration.
  The trajectory is integrated numerically.
 
  It reads the following input from stdin:
 
  - mu: RTBP mass parameter
  - t1: integration time
  - x[DIM]: initial condition, x=(l,L,g,G)
 
  The final point is written to stdout, in the form x=(l,L,g,G)
  
 */

int main( )
{
   double mu;
   double t1;			/* final time */
   double x[DIM];
   int status;

   // Input mass parameter, integration time and initial condition from stdin.
   if(scanf("%le %le %le %le %le %le", &mu, &t1, x, x+1, x+2, x+3) < 6)
   {
      perror("main: error reading input");
      exit(EXIT_FAILURE);
   }

   // Integrate trajectory numerically. 
   status=frtbp_del(mu,t1,x);
   if(status)
   {
      fprintf(stderr, "main: error integrating trajectory");
      exit(EXIT_FAILURE);
   }

   // Output final point to stdout.
   if(printf("%.15le %.15le %.15le %.15le\n", x[0], x[1], x[2], x[3])<0)
   {
      perror("main: error writting output");
      exit(EXIT_FAILURE);
   }
   exit(EXIT_SUCCESS);
}
