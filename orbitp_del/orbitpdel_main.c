// ======================================
// Orbit of Poincare map in Delaunay RTBP
// ======================================
// FILE:          $RCSfile: orbitpdel_main.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2011-03-01 13:48:05 $
//
// PURPOSE
// =======
// Consider the Restricted Three Body Problem in Delaunay variables.
// Let ${g=0} be a Poincare section, and $P(l,L,g,G)$ be the associated
// Poincare map.
// Given a point $p=(l,L,g,G)$, this function computes the orbit ${P^k(p)}$
// for $k=1,2,\dots,n$.
//
// OVERALL METHOD
// ==============
//
// 1. Input parameters from stdin:
// 
//    - mass parameter 
//    - number of desired points in the orbit
//
// 2. For each initial point $x$ do
//    2.1. Compute orbit of $x$ under the Poincare map: ${P^k(x)}$.
//    2.2. Output orbit of $x$ to stdout.

#include <stdio.h>
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE
#include <gsl/gsl_errno.h>	// gsl_set_error_handler_off
#include <rtbp.h>	// DIM
#include "orbitpdel.h"	// orbitp_del, print_orbit_del

int main( )
{
   double mu;
   double x[DIM];
   int npt;
   double *orbit;

   // Input mass parameter and number of points in the orbit from stdin.
   if(scanf("%le %d", &mu, &npt) < 2)
   {
      perror("main: error reading input");
      exit(EXIT_FAILURE);
   }

   orbit = (double *)malloc(npt*DIM*sizeof(double));

   // Stop GSL default error handler from aborting the program
   gsl_set_error_handler_off();

   while(scanf("%le %le %le %le", x, x+1, x+2, x+3)==4)
   {
      orbitp_del(mu,x,npt,orbit);
      print_orbit_del(npt,orbit);
   }
   exit(EXIT_SUCCESS);
}
