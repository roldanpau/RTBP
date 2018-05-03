/*! \file
    \brief Flow of the Reduced Restricted Three Body Problem: main prog.

    $Author: roldan $
    $Date: 2013-03-26 22:18:24 $
*/

#include <stdio.h>
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE
#include "frtbpred.h"

/**
  Flow of the Reduced Restricted Three Body Problem

  Compute the flow $\phi(s,x)$ of the reduced RTBP for a given time and
  initial condition. Reduced means that we identify $l$ with time and then
  we do not need to integrate this variable.
  The time $s$ may be positive or negative, allowing for forward or backward
  integration.
  The trajectory is integrated numerically.

  It reads the following input from stdin:
 
  - mu: RTBP mass parameter
  - s1: integration time
  - x[DIMRED]: Initial condition, 6 coordinates: (l, L, g, G, t, I).
 
  The final point is written to stdout, in the form 

  - x=(l,L,g,G,t,I)
  
 */
 
int main( )
{
   double mu;
   double s1;			/* final time */
   double x[DIMRED];
   int status;

   // Input mass parameter, integration time and initial condition from stdin.
   if(scanf("%le %le %le %le %le %le %le %le", &mu, &s1, x, x+1, x+2, x+3, x+4, x+5) < 8)
   {
      perror("main: error reading input");
      exit(EXIT_FAILURE);
   }

   // Integrate trajectory numerically. 
   status=frtbp_red_g(mu,s1,x);
   if(status)
   {
      fprintf(stderr, "main: error integrating trajectory");
      exit(EXIT_FAILURE);
   }

   // Output final point to stdout.
   if(printf("%.15le %.15le %.15le %.15le %.15le %.15le\n", x[0], x[1], x[2], x[3], x[4], x[5])<0)
   {
      perror("main: error writting output");
      exit(EXIT_FAILURE);
   }
   exit(EXIT_SUCCESS);
}
