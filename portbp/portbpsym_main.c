// ===========================================
// Find a symmetric periodic orbit of the RTBP
// ===========================================
// FILE:          $RCSfile: portbpsym_main.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2012-11-26 10:30:51 $
//
// PURPOSE:
// This program refines a trajectory of the RTBP which is close to a periodic
// orbit, until a true periodic orbit is obtained.
// We look for a periodic orbit in the given energy manifold $H=H_0$.
// Consider the Poincare section "sec" in the autonomous
// RTBP~\eqref{eq:RTBP}. We look for the periodic periodic orbit as a fixed
// point of the $k$-th iterate $P^{k}(x)$ of the Poincare first return map to
// this section.
// To refine the (approximate) fixed point, we use a Newton method.
// The program reads an approximate fixed point (x,px) from stdin. A Newton
// correction is applied to the initial condition until the refined initial
// condition (x,px) gives rise to a true fixed point of the 2D map
// $P^{k}(x)$.
// The refined fixed point is written to stdout.
// 
// Moreover, as a side product, the period of the periodic orbit is also
// written to stdout.
//
// NOTES:
//
// OVERALL METHOD:
//
// 1. Input parameters and approximate fixed point from stdin.
// 2. Refine trajectory of the RTBP which is close to a periodic orbit,
//    until a true periodic orbit is obtained.
// 3. Output true fixed point (gives rise to a true periodic orbit) and the
//    period of the periodic orbit.

#include <stdio.h>
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE
#include <string.h>     // strcmp
#include <gsl/gsl_errno.h>      // gsl_set_error_handler_off
#include <section.h>
#include <prtbp_2d.h>
#include "portbpsym.h"

int main( )
{
   double mu, H;
   section_t sec;
   int k;		// number of cuts with Poincare section
   double T2;		// half period of periodic orbit
   double x, x2[2];

   // auxiliary variables
   int status;
   char section_str[10];        // holds input string "SEC1", "SEC2" etc

   // Input mass parameter, Poincare section, energy value, number of cuts, 
   // initial condition from stdin.
   if(scanf("%le %s %le %d %le", &mu, section_str, &H, &k, &x) < 5)
   {
      perror("main: error reading input");
      exit(EXIT_FAILURE);
   }

   if (strcmp(section_str,"SEC1") == 0)
      sec = SEC1;
   else if (strcmp(section_str,"SEC2") == 0)
      sec = SEC2;
   else
   {
      perror("main: error reading section string");
      exit(EXIT_FAILURE);
   }

   // Stop GSL default error handler from aborting the program
   gsl_set_error_handler_off();

   // Refine trajectory of the RTBP which is close to a periodic orbit, until
   // a true periodic orbit is obtained.
   status = portbpsym(mu,sec,H,k,&x);
   if(status)
   {
      fprintf(stderr, \
	    "main: unable to find periodic orbit up to desired accuracy\n");
      // this is not really an error
      // exit(EXIT_FAILURE);
   }
   // Compute half period T2 by calling prtbp_2d
   x2[0]=x;
   x2[1]=0;	// px
   if(prtbp_2d(mu,sec,H,k,x2,&T2))
   {
      fprintf(stderr, "dist: error computing 2D poincare map\n");
      exit(EXIT_FAILURE);
   }
   // Output root and half period to stdout.
   if(printf("%.15le %.15le\n", x, T2)<0)
   {
      perror("main: error writting output");
      exit(EXIT_FAILURE);
   }
   exit(EXIT_SUCCESS);
}
