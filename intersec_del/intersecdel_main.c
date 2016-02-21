// ===================================
// Intersection of Invariant Manifolds
// ===================================
// FILE:          $RCSfile: intersecdel_main.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2012-06-06 11:07:40 $
//
// PURPOSE
// =======
// Consider the Restricted Three Body Problem in Delaunay variables. 
// Assume that energy is fixed: $H(l,L,g,G)=\bar H$.
// Let SEC2= \{l=\pi\}$ be the Poincare section, and $\sixmap$ be the
// ``iterated'' 2D Poincare map.
// Let $p$ be a fixed point of $\sixmap$ associated to the 3:1 resonant
// periodic orbit of the flow, and let $W^u(p)$ be the associated unstable
// manifold.
// This function computes the ``primary'' intersection point $z=(l,G)$ with
// the axis line defined by
//    \[ g = l. \]
// Due to the symmetry of the manifolds in this section, it is enough to look
// for the first intersection of the unstable manifold $W^u(p)$ with the
// reversibility line (the $G$ axis of the poincare section).
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
//    - eccentricity "e"
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
//    - eccentricity e
//    - energy level H
//    - point p_u,
//    - integration time t_u to reach the intersection point z, 
//    - intersection point z = P(p_u).

#include <stdio.h>
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE
#include "intersecdel.h"

void print_pt(double z[2])
{
      if(printf("%.15le %.15le ", z[0], z[1])<0)
      {
         perror("main: error writting output");
         exit(EXIT_FAILURE);
      }
}

int main( )
{
   double mu, e, H;

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

   double p_u[2];	// point in the unstable segment
   double z[2];		// homoclinic point

   double t;		// integration time to reach z from p_u/p_s

   // auxiliary vars
   int status;

   // 1. Input parameters from stdin.
   if(scanf("%le %d %le", &mu, &stable, &l) < 3)
   {
      perror("main: error reading input");
      exit(EXIT_FAILURE);
   }

   while(scanf("%le %le %le %le %le %le %le %d %le %le", 
	    &e, &H, p, p+1, v, v+1, &lambda, &n, &h1, &h2) == 10)
   {
      // 2. Find a root of the distance function, i.e. an intersection point of
      // the manifolds
      if(!stable)
      {
	 status = intersec_del_unst(mu, H, p, v, lambda, n, h1, h2, l, &h, p_u, &t, z);
	 if(status)
	 {
	    fprintf(stderr, "main: error computing intersection point\n");
	    exit(EXIT_FAILURE);
	 }
      }
      else
      {
	 status = intersec_del_st(mu, H, p, v, lambda, n, h1, h2, l, &h, p_u, &t, z);
	 if(status)
	 {
	    fprintf(stderr, "main: error computing intersection point\n");
	    exit(EXIT_FAILURE);
	 }
      }

      // 3. Output the following data to stdout:
      //    - eccentricity e
      //    - energy level H
      //    - point p_u/p_s
      //    - integration time t to reach the intersection point z, 
      //    - intersection point z = P(p_u).
      printf("%e ", e);
      printf("%.15le ", H);
      print_pt(p_u);
      printf("%.15le ", t);
      print_pt(z);
      printf("\n");
   }
   exit(EXIT_SUCCESS);
}
