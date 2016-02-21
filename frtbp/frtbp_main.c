// FILE:          frtbp_main.c
// TITLE:         Flow of the Restricted Three Body Problem
// AUTHOR:        Pau Roldan
//                All code is my own except where credited to others.
// DATE:          September 30, 2009
//
// PURPOSE:
// Compute the flow $\phi(t,x)$ of the RTBP for a given time (positive or
// negative) and initial condition.
// It reads the following input from stdin:
//
// - mu: RTBP mass parameter
// - t1: integration time
// - x[DIM]: initial condition, x=(X,Y,PX,PY)
//
// The trajectory is integrated numerically. The final point is written to
// stdout, in the form x=(X,Y,PX,PY).
//
// NOTES:
// Integration method: Taylor method (provided by Angel Jorba).
//
// OVERALL METHOD:
// The list of general tasks is:
// 1. Input integration time and initial condition from stdin.
// 2. Integrate trajectory numerically. 
// 3. Otput final point to stdout.
//
// FUNCTIONS:

#include <stdio.h>
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE
#include <math.h>
#include "frtbp.h"

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
   status=frtbp(mu,t1,x);
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
