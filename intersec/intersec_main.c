// ===================================
// Intersection of Invariant Manifolds
// ===================================
// FILE:          $RCSfile: intersec_main.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2013-02-20 10:00:58 $
//
// PURPOSE
// =======
// Consider the Restricted Three Body Problem in euclidean variables. 
// Assume that energy is fixed: $H(x,y,p_x,p_y)=\bar H$.
// Let $\Sigma'=\{y=0, p_y<0\}$ be the Poincare section, and $\sixmap$ be the
// ``iterated'' 2D Poincare map.
// Let $p_3$ be a fixed point of $\sixmap$ associated to the 7:1 resonant
// periodic orbit of the flow.
// Let $W^u(p_3)$ be the associated unstable manifold.
// This function computes the ``primary'' intersection point $z=(x,h)$ with
// the axis line defined by
//    \[ p_x = l. \]
//
// OVERALL METHOD
// ==============
//
// 1. Input parameters from stdin:
// 
//    - mass parameter 
//    - stable flag
//    - axis line "l$
//
// For each energy value, do:
//
// 2.1 Input data:
//    - energy value "H"
//    - fixed point "p"
//    - linear unstable direction "v"
//    - unstable eigenvalue "lambda"
//    - number of iterations "n" in the unstable dir by the Poincare map to
//       reach the intersection point
//    - "h1, h2" small increment in the direction of v (initial guess for
//    bisection).
//
// 2.2. Find a root of the distance function, i.e. an intersection point of the
// manifold with the axis line.
//
// 2.3. Output the following data to stdout:
//    - energy level H
//    - point p_u,
//    - integration time t_u to reach the intersection point z, 
//    - intersection point z = P(p_u).

#include <stdio.h>
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE
#include <rtbp.h>	// DIM
#include <utils_module.h>	// dblprint

#include "intersec.h"

int main( )
{
   double mu, H;

   double p[2];		// fixed point
   double v[2];		// eigenvector
   double lambda;	// eigenvalue

   // number of desired iterations by the Poincare map
   int n;		// in the unstable dir

   // Small increment in the direction of v (initial guess for Newton).
   // This is furnished by program "approxint".
   double h1, h2;

   // "stable" flag specifies wheather we want to compute the unstable (=0)
   // or stable (=1) manifold
   int stable;

   double l;		// axis line p_x = l

   double h;		// root of distance function

   double p_u[DIM];	// point in the unstable segment
   double z[DIM];	// homoclinic point

   double t;		// integration time to reach z from p_u/p_s

   // auxiliary vars
   int status;

   // 1. Input parameters from stdin.
   if(scanf("%le %d %le", &mu, &stable, &l) < 3)
   {
      perror("main: error reading input");
      exit(EXIT_FAILURE);
   }

   while(scanf("%le %le %le %le %le %le %d %le %le", 
	    &H, p, p+1, v, v+1, &lambda, &n, &h1, &h2) == 9)
   {
      // 2. Find a root of the distance function, i.e. an intersection point of
      // the manifolds
      if(!stable)
      {
	 status = intersec_unst(mu, H, p, v, lambda, n, h1, h2, l, &h, p_u, &t, z);
	 if(status)
	 {
	    fprintf(stderr, "main: error computing intersection point\n");
	    exit(EXIT_FAILURE);
	 }
      }
      else
      {
	 status = intersec_st(mu, H, p, v, lambda, n, h1, h2, l, &h, p_u, &t, z);
	 if(status)
	 {
	    fprintf(stderr, "main: error computing intersection point\n");
	    exit(EXIT_FAILURE);
	 }
      }

      // 3. Output the following data to stdout:
      //    - energy level H
      //    - point p_u/p_s
      //    - integration time t to reach the intersection point z, 
      //    - intersection point z = P(p_u).
      printf("%.15le ", H);
      dblprint(p_u,DIM);
      printf("%.15le ", t);
      dblprint(z,DIM);
      printf("\n");
   }
   exit(EXIT_SUCCESS);
}
