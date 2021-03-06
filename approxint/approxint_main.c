/*! \file
    \brief Approximate Intersection of Invariant Manifolds: main prog
    \author Pau Roldan
*/

#include <stdio.h>
#include <stdlib.h>		// EXIT_SUCCESS, EXIT_FAILURE
#include <section.h>	// branch_t
#include <prtbp.h>		// SEC2
#include <errmfld.h>	// h_opt
#include "approxint.h"	// approxint_unst, approxint_st

/** 
   Approximate Intersection of Invariant Manifolds: main prog

   Consider the Restricted Three Body Problem. Assume that energy is fixed:
   \f$ H(x,y,p_x,p_y)=\bar H \f$.
   Let SEC2={y=0, p_y<0} be a Poincare section, and $P(x,p_x)$ be 
   the associated 2D Poincare map.
   Let $p(x,p_x)$ be a fixed point of the Poincare map associated to the 1:3
   resonant periodic orbit of the flow.
   Let $W^u(p),W^s(p)$ be the unstable,stable manifold of $p$.
   This function computes an approximation to the "first" intersection of the
   manifold with the $x$ axis, or more generally with any horizontal line \f$
   p_x=a. \f$
  
   OVERALL METHOD
  
   1. Input parameters from stdin:
   
      - mass parameter 
      - number of cuts "k" with Poincare section
      - "stable" flag (unstable=0, stable=1)
      - "branch" flag (left=0, right=1)
      - axis line "a"
  
   2. For each input line, do
   
      2.1. Input parameters from stdin:
         - energy value "H"
         - fixed point "p"
         - linear unstable,stable direction "v"
         - unstable,stable eigenvalue "lambda"
  
      Estimate error commited in the linear approximation of the manifold
  
      2.2. Find an approximate intersection point of the manifolds, in the
      form of an interval $u_i=(h_1,h_2)$ containig the root p_u = p+h_u v_u,
      where \f$h_u\in (h_1,h_2)\f$.
  
      2.3. Output the following line to stdout: 
         H, iter, h_1, h_2, z.

  \pre Initialize the aproximation to the first homoclinic point by hand!
 */

int main( )
{
   double mu, H;
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

   // "branch" flag specifies wheather we want to compute the left branch (=0)
   // or right branch (=1) of the manifold
   int branch;

   // approximate homoclinic point
   double z[2];

   // H=-1.3594 (stable manifold)
   //double z[2] = {-0.0438623, 0};
   // H=-1.3594 (unstable manifold)
   //double z[2] = {-0.0669547, 0};
   // H=-1.343403e+00
   //double z[2] = {-0.0607277, 0};
   // H=-1.400403e+00
   //double z[2] = {-8.464640048727097e-02, 0};
   // H=-1.7194 (stable manifold)
   //double z[2] = {-0.3647, 0};
   // H=-1.7194 (unstable manifold)
   //double z[2] = {-0.4066, 0};

   double a;	// horizontal axis line $p_x=a$

   // auxiliary vars
   int status;
   branch_t br;

   // 1. Input parameters from stdin.
   if(scanf("%le %d %d %d %le", &mu, &k, &stable, &branch, &a) < 5)
   {
      perror("main: error reading input");
      exit(EXIT_FAILURE);
   }

   br = (branch==0 ? LEFT : RIGHT);

   // Initialize approximation to homoclinic point by hand
   if(!stable)
   {
      if(br==RIGHT)
      {
	 z[0] = -0.4066;
	 z[1] = 0;
      }
      else //br==LEFT
      {
	 z[0] = -0.5902;
	 z[1] = 0;
      }
   }
   else
   {
      if(br==RIGHT)
      {
	 z[0] = -0.3647;
	 z[1] = 0;
      }
      else //br==LEFT
      {
	 z[0] = -0.5626;
	 z[1] = 0;
      }
   }

   // For each energy level H in the range, do
   while(scanf("%le %le %le %le %le %le", 
	    &H, p, p+1, v, v+1, &lambda)==6)
   {
      fprintf(stderr, "\nH: %e\n", H);

	  h = 1.e-2;	// upper bound on linear displacement

      // By default we work with the RIGHT branch of the manifolds.
      // To work with the LEFT branch of the manifolds instead,
      // we just take the negative of the displacement h in the linear
      // approximation.
      if(br==LEFT) h=-h;

      // Compute optimal displacement h such that the estimate error commited
      // in the linear approximation of the manifold is smallest.
      h = h_opt(mu,SEC2,H,k,p,v,lambda,stable,h);

      // 2. Find an approximate intersection point of the unstable manifold and
      // the $x$ axis (approximate manifold by segments, and check if they
      // cross the $x$ axis).
      if(!stable)
      {
	 status = approxint_unst(mu, H, k, p, v, lambda, h, a, &iter, &h_1, 
	       &h_2, z);
	 if(status)
	 {
	    fprintf(stderr, 
		  "H=%e: couldn't find approx. intersection point\n", H);
	    //exit(EXIT_FAILURE);
	 }
	 else
	 {
		  // 3. Output the following data to stdout:
		  //    H, h_1, h_2, z
		  printf("%.15e %d %.15e %.15e %.15e %.15e\n", H, iter, h_1, h_2, z[0], z[1]);
		  fflush(NULL);
	 }
      }
      else
      {
	 status = approxint_st(mu, H, k, p, v, lambda, h, a, &iter, &h_1, 
	       &h_2, z);
	 if(status)
	 {
	    fprintf(stderr, 
		  "H=%e: couldn't find approx. intersection point\n", H);
	    exit(EXIT_FAILURE);
	 }
	 else
	 {
		  // 3. Output the following data to stdout:
		  //    H, h_1, h_2, z
		  printf("%.15e %d %.15e %.15e %.15e %.15e\n", H, iter, h_1, h_2, z[0], z[1]);
		  fflush(NULL);
	 }
      }
      

   }
   exit(EXIT_SUCCESS);
}
