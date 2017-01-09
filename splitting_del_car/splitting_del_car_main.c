// ======================================
// Splitting Angle of Invariant Manifolds
// ======================================
// FILE:          $RCSfile: splitting_del_car_main.c,v $
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
//    - Poincare section "sec" where the homoclinic point lies
//    - linear unstable direction "v_u" in section SEC1
//    - linear unstable direction "v_u" in section SEC2
//    - number of iterations "n" in the unstable dir by the Poincare map to
//       reach the fold
//    - points in the unstable segment
//
// 2. Find splitting angle
//
// 3. Output the following data to stdout:

#include <stdio.h>
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE
#include <string.h>	// strcmp
#include <math.h>	// M_PI

#include <section.h>            // section_t
#include <sec2sec1_module.h>    // sec2sec1

#include "splitting_del_car.h"

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
   section_t sec;

   double v_u_SEC1[2];	// unstable vector on section SEC1
   double v_u_SEC2[2];	// unstable vector on section SEC2

   // number of desired iterations by the Poincare map
   int n;		// in the unstable dir

   double p_u[2];	// points in the unstable segment
   double angle;	// splitting angle

   // auxiliary vars
   double v_u[2];           // unstable vector
   int status;
   char section_str[10];    // holds input string "SEC1", "SEC2" etc
   double ti;

   // 1. Input parameters from stdin.
   if(scanf("%le", &mu) < 1)
   {
      perror("main: error reading input");
      exit(EXIT_FAILURE);
   }

   while(scanf("%le %s %le %le %le %le %d %le %le", 
	    &H, section_str, v_u_SEC1, v_u_SEC1+1, v_u_SEC2, v_u_SEC2+1, 
        &n, p_u, p_u+1) == 9)
   {
       if (strcmp(section_str,"SEC1") == 0)
          sec = SEC1;
       else if (strcmp(section_str,"SEC2") == 0)
          sec = SEC2;
       else
       {
          perror("main: error reading section string");
          exit(EXIT_FAILURE);
       }

       if(sec==SEC1)
       {
           // take p_u to section SEC1
           sec2sec1(mu,H,p_u,&ti);
           v_u[0]=v_u_SEC1[0];
           v_u[1]=v_u_SEC1[1];
       }
       else
       {
           v_u[0]=v_u_SEC2[0];
           v_u[1]=v_u_SEC2[1];
       }
       
      // 2. Find splitting angle
      status = splitting_del_car_unst(mu, H, sec, v_u, n, p_u, &angle);
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
