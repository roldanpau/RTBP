// ===============================================================
// Restricted Three Body Problem equations in Delaunay coordinates
// ===============================================================
// FILE:          rtbpdel_main.c
// AUTHOR:        Pau Roldan
//                All code is my own except where credited to others.
// DATE:          September 30, 2009
//
// PURPOSE:
//    The program computes the vectorfield of the RTBP at a point named "x".
//    Specifically, we use the equations of motion of the planar circular
//    RTBP in Delaunay coordinates (l,L,g,G), in the rotating (synodic)
//    coordinate system.
//    This system is Hamiltonian, with Hamiltonian function
//       H(l,L,g,G) = -1/(2L^2) - G + R,
//    where R is the perturbation relative to the two body problem.
//    See Szebehelly, page 364.
//    See Marcel's notes "Inner and outer dynamics".
//    See also my hand-written notes (yellow notebook).
//
// It reads the following input from stdin:
//
// - mu: RTBP mass parameter
// - x[DIM]: point in phase space, 4 coordinates: (l,L,g,G)
//
// The vectorfield at (t,x) is written to stdout, 4 coordinates: d/dt(l, L,
// g, G).
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
#include "rtbpdel.h"

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
   status=rtbp_del(0.0,x,y,&mu);
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
