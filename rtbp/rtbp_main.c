// FILE:          frtbp_main.c
// TITLE:         Flow of the Restricted Three Body Problem
// AUTHOR:        Pau Roldan
//                All code is my own except where credited to others.
// DATE:          September 30, 2009
//
// PURPOSE:
//    Compute the vectorfield of the RTBP at a point named "x".
//    Specifically, we use the equations of motion of the spatial circular
//    RTBP, in the rotating (synodic) coordinate system, using adimensional
//    quantities, and Hamiltonian formulation (position-momentum).
//    See Szebehelly, page 349. 
//    See also the notes on normal forms by Maciej Capinski.
// It reads the following input from stdin:
//
// - mu: RTBP mass parameter
// - x[DIM]: point in phase space, 4 coordinates: (X, Y, P_X, P_Y).
//
// The vectorfield at (t,x) is written to stdout, 4 coordinates: d/dt(X, Y,
// P_X, P_Y).
//
// NOTES:
//
// OVERALL METHOD:
// The list of general tasks is:
// 1. Input mass parameter and point from stdin.
// 2. Compute vectorfield
// 3. Otput vectorfield to stdout.
//
// FUNCTIONS:

#include <stdio.h>
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE
#include <math.h>
#include "rtbp.h"

int main( )
{
   double mu;
   double x[DIM];
   double y[DIM];
   int status;

   // Input mass parameter, and phase space point from stdin.
   if(scanf("%le %le %le %le %le", &mu, x, x+1, x+2, x+3) < 5)
   {
      perror("main: error reading input");
      exit(EXIT_FAILURE);
   }

   // Compute vectorfield
   status=rtbp(0.0,x,y,&mu);
   if(status)
   {
      fprintf(stderr, "main: error computing vectorfield");
      exit(EXIT_FAILURE);
   }

   // Output vectorfield to stdout.
   if(printf("%.15le %.15le %.15le %.15le\n", y[0], y[1], y[2], y[3])<0)
   {
      perror("main: error writting output");
      exit(EXIT_FAILURE);
   }
   exit(EXIT_SUCCESS);
}
