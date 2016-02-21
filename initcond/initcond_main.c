// =============================================
// initial condition of p/q resonant p.o. in 2BP
// =============================================
// FILE:          $RCSfile: initcond_main.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2011-04-12 10:11:03 $
//
// PURPOSE:
// Consider the Asteroid - Sun/Jupiter two body problem (2BP).
// Given an energy value "h", this program computes an initial condition of
// Asteroid (x,y,p_x,p_y) corresponding to r=p/q ressonant periodic orbit in
// the 2BP.
//
// We impose that initial condition is at the perihelion. This implies that
// position of Asteroid is A=(x,0) and velocity (and thus momentum) is
// v=(0,v_y). Therefore, we output as initial condition only the components
// (x,v_y).
//
// NOTES:
//
// OVERALL METHOD:
//
// 1. Input parameters (energy level H, resonance r) from stdin.
// 2. Compute initial condition.
// 3. Output initial condition.

#include <stdio.h>
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE
#include "initcond.h"

int main( )
{
   double H, r, x, py;
   int status;

   // Input energy value, resonance from stdin.
   if(scanf("%le %le", &H, &r) < 2)
   {
      perror("main: error reading input");
      exit(EXIT_FAILURE);
   }

   // Compute initial condition.
   status = initcond(H,r,&x,&py);
   if(status)
   {
      fprintf(stderr, \
	    "main: unable to compute initial condition\n");
      exit(EXIT_FAILURE);
   }
   // Output initial condition to stdout.
   if(printf("%.15le %.15le\n", x, py)<0)
   {
      perror("main: error writting output");
      exit(EXIT_FAILURE);
   }
   exit(EXIT_SUCCESS);
}
