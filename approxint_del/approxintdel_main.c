// ===============================================
// Approximate Intersection of Invariant Manifolds
// ===============================================
// FILE:          $RCSfile: approxintdel_main.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2012-06-06 10:15:20 $
//
// PURPOSE
// =======
// Consider the Restricted Three Body Problem. Assume that energy is fixed:
// $H(l,L,g,G)=\bar H$.
// Let SEC2: {l=\pi} be a Poincare section, and $P(g,G)$ be the associated 2D
// Poincare map.
// Let $p(g,G)$ be a fixed point of the Poincare map associated to the 3:1
// resonant periodic orbit of the flow.
// Let $W^u(p),W^s(p)$ be the unstable,stable manifold of $p$.
// This function computes an approximation to the "first" intersection of the
// manifold with the $G$ axis.
//
// OVERALL METHOD
// ==============
//
// 1. Input parameters from stdin:
// 
//    - mass parameter 
//    - number of cuts "k" with Poincare section
//    - stable flag
//    - axis line "a"
//
// 2. For each input line, do
//    2.1. Input parameters from stdin:
//       - eccentricity "e"
//       - energy value "H"
//       - fixed point "p"
//       - linear unstable,stable direction "v"
//       - unstable,stable eigenvalue "lambda"
//
//    Estimate error commited in the linear approximation of the manifold
//
//    2.2. Find an approximate intersection point of the manifolds, in the
//    form of an interval $u_i=(h_1,h_2)$ containig the root p_u = p+h_u v_u,
//    where $h_u\in (h_1,h_2)$.
//
//    2.3. Output the following line to stdout: 
//       e, H, iter, h_1, h_2, z.

#include <stdio.h>
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE
#include <errmflddel.h>	// h_opt_del

#include "approxintdel.h"	// approxint_del_unst, approxint_del_st

int main( )
{
   double mu, e, H;
   int k;

   double p[2];		// fixed point
   double v[2];		// unstable/stable vector
   double lambda;	// unstable/stable eigenvalue

   // small increment in the direction of v
   double h;

   // Number of iterations of poinc map to reach homoclinic point.
   int iter;

   // Displacements in the unst direction straddling the root h.
   double h_1, h_2;

   // "stable" flag specifies wheather we want to compute the unstable (=0)
   // or stable (=1) manifold
   int stable;

   // Approximate homoclinic point, obtained by visual inspection of
   // invariant manifolds.  e=0.84, lower separatrix
   double z[2] = {0.0, 0.34528};

   double a;	// vertical axis line $G=a$

   // auxiliary vars
   int status;

   // 1. Input parameters from stdin.
   if(scanf("%le %d %d %le", &mu, &k, &stable, &a) < 4)
   {
      perror("main: error reading input");
      exit(EXIT_FAILURE);
   }

   // For each eccentricity level "e" in the range, do
   while(scanf("%le %le %le %le %le %le %le", 
	    &e, &H, p, p+1, v, v+1, &lambda)==7)
   {
      fprintf(stderr, "\ne: %e\n", e);

      // Compute optimal displacement h such that the estimate error commited
      // in the linear approximation of the manifold is smallest.
      h = h_opt_del(mu,H,k,p,v,lambda,stable);

      // 2. Find an approximate intersection point of the unstable manifold and
      // the $G$ axis (approximate manifold by segments, and check if they
      // cross the $G$ axis).
      if(!stable)
      {
	 status = approxint_del_unst(mu, H, k, p, v, lambda, h, a, &iter,
	       &h_1, &h_2, z);
	 if(status)
	 {
	    fprintf(stderr, 
		  "e=%e: couldn't find approx. intersection point\n", e);
	    // 3. Output the following data to stdout:
	    //    e, H, iter, h_1, h_2, z
	    printf("%e %.15e %d %.15e %.15e %.15e %.15e\n", e, H, 0, 0.0, 0.0, 0.0, 0.0);
	    fflush(NULL);	// flush all open output streams
	    continue;
	 }
      }
      else
      {
	 status = approxint_del_st(mu, H, k, p, v, lambda, h, a, &iter, 
	       &h_1, &h_2, z);
	 if(status)
	 {
	    fprintf(stderr, 
		  "e=%e: couldn't find approx. intersection point\n", e);
	    // 3. Output the following data to stdout:
	    //    e, H, iter, h_1, h_2, z
	    printf("%e %.15e %d %.15e %.15e %.15e %.15e\n", e, H, 0, 0.0, 0.0, 0.0, 0.0);
	    fflush(NULL);	// flush all open output streams
	    continue;
	 }
      }
      

      // 3. Output the following data to stdout:
      //    e, H, iter, h_1, h_2, z
      printf("%e %.15e %d %.15e %.15e %.15e %.15e\n", e, H, iter, h_1, h_2, z[0], z[1]);
      fflush(NULL);	// flush all open output streams
   }
   exit(EXIT_SUCCESS);
}
