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
//    - section type "sec": sec={SEC1,SEC2}
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
//    - Point in local unstable manifold of the appropriate pendulum
//    (Delaunay coords).
 

#include <stdio.h>
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE
#include <string.h>	// strcmp
#include <rtbp.h>	// DIM

#include <section.h>
#include "intersec_del_car.h"

#include <utils_module.h>	// dblprint

int main( )
{
   double mu, H;
   section_t sec;   // Poincare section

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

   double z_car[DIM];		// homoclinic point in Cartesian
   double z_del[DIM];		// homoclinic point in Delaunay
   double z_u[DIM];		// point in local unstable manifold

   double t;		// integration time to reach z from p_u/p_s

   // auxiliary vars
   char section_str[10];    // holds input string "SEC1", "SEC2" etc
   int status;

   // 1. Input parameters from stdin.
   if(scanf("%le %s %d %le", &mu, section_str, &stable, &l) < 4)
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

   while(scanf("%le %le %le %le %le %le %d %le %le", 
	    &H, p, p+1, v, v+1, &lambda, &n, &h1, &h2) == 9)
   {
      // 2. Find a root of the distance function, i.e. an intersection point of
      // the manifolds
      if(!stable)
      {
          status = intersec_del_car_unst(mu, sec, H, p, v, lambda, n, 
                  h1, h2, l, &h, p_u, &t, z_del, z_car, z_u);
          if(status)
          {
              fprintf(stderr, "main: error computing intersection point\n");
              exit(EXIT_FAILURE);
          }
      }
      else
      {
	 status = intersec_del_car_st(mu, H, p, v, lambda, n, h1, h2, l, &h, p_u, &t, z);
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
      //    - point in local unstable manifold z_u
      printf("%.15le ", H);
      dblprint(p_u,2);
      printf("%.15le ", t);
      dblprint(z_del,DIM);
      dblprint(z_car,DIM);
      dblprint(z_u,DIM);
      printf("\n");
      fflush(NULL);
   }
   exit(EXIT_SUCCESS);
}
