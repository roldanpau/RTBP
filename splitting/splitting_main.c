// ======================================
// Splitting Angle of Invariant Manifolds
// ======================================
// FILE:          $RCSfile: splitting_main.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2013-03-11 11:36:14 $
//
// PURPOSE
// =======
// Consider the Restricted Three Body Problem. Assume that energy is fixed:
// $H(x,y,p_x,p_y)=\bar H$.
// Let ${y=0} be a Poincare section, and $P(x,p_x)$ be the associated 2D
// Poincare map.
// Let $p(x,p_x)$ be a fixed point of the Poincare map associated to a
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
//    - energy value "H"
//    - linear unstable direction "v_u"
//    - number of iterations "n" in the unstable dir by the Poincare map to
//       reach the fold
//    - points in the unstable segment
//
// 2. Find splitting angle
//
// 3. Output the following data to stdout:

#include <stdio.h>
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE
#include <math.h>	// M_PI
#include "splitting.h"

void print_pt(double z[2])
{
      if(printf("% .15le % .15le\n", z[0], z[1])<0)
      {
         perror("main: error writting output");
         exit(EXIT_FAILURE);
      }
}

int main( )
{
   double mu, H;

   double v_u[2];	// unstable vector

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

   while(scanf("%le %le %le %d %le %le", 
	    &H, v_u, v_u+1, &n, p_u, p_u+1) == 6)
   {
      // 2. Find splitting angle
      status = splitting_st(mu, H, v_u, n, p_u, &angle);
      if(status)
      {
	 fprintf(stderr, "main: error computing splitting angle");
	 exit(EXIT_FAILURE);
      }

      // 3. Output the following data to stdout:
      //    - splitting angle (in radians)
      printf("%.15le %.15le\n", H, angle);
   }
   exit(EXIT_SUCCESS);
}
