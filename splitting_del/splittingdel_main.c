// ======================================
// Splitting Angle of Invariant Manifolds
// ======================================
// FILE:          $RCSfile: splittingdel_main.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2012-06-06 13:07:51 $
//
// PURPOSE
// =======
// Consider the Restricted Three Body Problem. Assume that energy is fixed:
// $H(l,L,g,G)=\bar H$.
// Let SEC2: {l=pi} be a Poincare section, and $P(g,G)$ be the associated 2D
// Poincare map.
// Let $p(g,G)$ be a fixed point of the Poincare map associated to a
// periodic orbit of the flow.
// Let $W^u(p), W^s(p)$ be the unstable/stable manifold of $p$.
// Let $z$ be the intersection point $z\in W^u(p) \cap W^s(p)$ that is closer
// to the "fold" of the manifolds. 
// This program computes the splitting angle between the manifolds at the
// intersection point $z$.
//
// OVERALL METHOD
// ==============
//
// 1. Input parameters from stdin:
// 
//    - mass parameter 
//    - eccentricity "e"
//    - energy value "H"
//    - linear unstable direction "v_u"
//    - number of iterations "n" in the unstable dir by the Poincare map to
//       reach the fold
//    - points in the unstable segment
//
// 2. Find splitting angle
//
// 3. Output the following data to stdout:
//    - tangent vector to the manifold w[2]
//    - splitting angle (in radians)

#include <stdio.h>
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE
#include <math.h>	// M_PI

#include "splittingdel.h"

int main( )
{
   double mu, e, H;

   double v_u[2];	// unstable vector
   double w_u[2];	// tangent vector to the manifold at hom.pt

   // number of desired iterations by the Poincare map
   int n;		// in the unstable dir

   double p_u[2];	// points in the unstable segment
   double angle;	// splitting angle

   // auxiliary vars
   int status;

   // 1. Input parameters from stdin.
   if(scanf("%le", &mu) < 1)
   {
      perror("main: error reading input");
      exit(EXIT_FAILURE);
   }

   while(scanf("%le %le %le %le %d %le %le", 
	    &e, &H, v_u, v_u+1, &n, p_u, p_u+1) == 7)
   {
      // 2. Find splitting angle
      status = splitting_del(mu, H, v_u, n, p_u, w_u, &angle);
      if(status)
      {
	 fprintf(stderr, "main: error computing splitting angle");
	 exit(EXIT_FAILURE);
      }

      // 3. Output the following data to stdout:
      //    - splitting angle (in radians)
      printf("%e %.15le %.15le %.15le %.15le\n", e, H, w_u[0], w_u[1], angle);
   }
   exit(EXIT_SUCCESS);
}
