// ======================================
// Splitting Angle of Invariant Manifolds
// ======================================
// FILE:          $RCSfile: axis_main.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2013-01-25 12:09:03 $
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
// 2. For each input line, do:
//
//    Input data:
//    - energy value "H"
//    - x_u: intersection of unstable manifold with the line $p_y=1.e-5$ (x
//    coordinate).
//    - x_s: intersection of stable manifold with the line $p_y=1.e-5$ (x
//    coordinate).
//
// 2. Find splitting angle using numerical differentiation, plus Richardson
// extrapolation.
//
// 3. Output the following data to stdout:
//    - H
//    - splitting angle

#include <stdio.h>
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE
#include <math.h>	// M_PI

int main( )
{
   const double h = 1.e-5;	// step size for central differences

   double mu, H;

   double xu_p1, xs_p1;	// intersection of manifold with line p_y=h
   double xu_m1, xs_m1;	// intersection of manifold with line p_y=-h
   double xu_p2, xs_p2;	// intersection of manifold with line p_y=h
   double xu_m2, xs_m2;	// intersection of manifold with line p_y=-h
   double extra;	// extrapolation
   double angle;	// splitting angle

   // auxiliary variables
   double d1,d2;	// central differences

   // 1. Input parameters from stdin.
   if(scanf("%le", &mu) < 1)
   {
      perror("main: error reading input");
      exit(EXIT_FAILURE);
   }

   while(scanf("%le %le %le %le %le %le %le %le %le", 
	    &H, &xu_p1, &xs_p1, &xu_m1, &xs_m1, 
	    &xu_p2, &xs_p2, &xu_m2, &xs_m2) == 9)
   /*
   while(scanf("%le %le %le %le %le", 
	    &H, &xu_p1, &xs_p1, &xu_m1, &xs_m1) == 5)
	    */
   {
      // 2. Find splitting angle
      d1 = ((xu_p1-xs_p1)-(xu_m1-xs_m1))/(2*h);
      d2 = ((xu_p2-xs_p2)-(xu_m2-xs_m2))/(4*h);
      extra = (4*d1-d2)/3.0;
      angle = atan(extra);
      //angle = atan(d1);

      // 3. Output the following data to stdout:
      //    - splitting angle (in radians)
      printf("%e %.15le %.15le %.15le %.15le\n", H, angle,d1,d2,extra);
   }
   exit(EXIT_SUCCESS);
}
