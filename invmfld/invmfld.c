/*! \file invmfld.c
    \brief Invariant Manifolds of fixed point of Poincare map in RTBP
    
    Consider the Restricted Three Body Problem. Assume that energy is fixed:
    \f$H(x,y,p_x,p_y)=\bar H\f$.
    Let "sec" be a Poincare section, and $P(x,p_x)$ be the associated 2D
    Poincare map (consisting of $k$ cuts with the Poincare section).
    Let $p(x,p_x)$ be a fixed point of the Poincare map associated to a
    periodic orbit of the flow.
    This function computes the unstable/stable manifolds $W^u(p), W^s(p)$ of
    the point $p$. A manifold is represented numerically as a sets of points.

    $Author: roldan $
    $Date: 2013-03-11 11:34:06 $
*/

#include <stdio.h>
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE
#include <string.h>	// strcmp
#include <math.h>	// sqrt
#include <lift.h>
#include <utils_module.h>	// dblprint
#include <prtbp_nl_2d_module.h>	// prtbp_nl_2d, prtbp_nl_2d_inv
#include <prtbp_nl.h>	// prtbp_nl, prtbp_nl_inv
#include <errmfld.h>
#include "disc.h"	// disc

/// Number of points in discretization of linear segment
const int NPOINTS = 100; 

/**
  Main program.

  Input params (stdin): 
  
//    - mass parameter 
//    - Poincare section "sec"
//    - energy value "H"
//    - number of cuts "k" with Poincare section
//    - fixed point "p"
//    - linear unstable (or stable) direction "v"
//    - unstable (or stable) eigenvalue "lambda"
//    - number of desired iterations by the Poincare map "n"
//    - "stable": flag specifying whether to compute the unstable manifold
//    (stable==0) or the stable manifold (stable==1).
//    - h: small increment in the direction of v
//    - ifp: Index specifying which iterate of fixed point we are interested
//    in, e.g. p_2 -> ifp==2

  Output params (stdout): sequence of points approximating the manifold.

  \remark
  If the flag "stable" specifies the unstable manifold (0), we iterate the
  forward Poincare map $P$.

  \remark
  If the flag "stable" specifies the stable manifold (1), we iterate the
  backward Poincare map $P^{-1}$.
*/

// NOTE: instead of passing h as a parameter, it would be better to call
// h_opt().

int main( )
{
   double mu, H;
   section_t sec;
   int k;		// number of iterates of Poincare map

   double p[2];		// fixed point
   double v[2];		// stable/unstable direction
   double lambda;	// stable/unstable eigenvalue
   int n;		// number of desired iterations of linear segment

   // "stable" flag specifies wheather we want to compute the unstable (=0)
   // or stable (=1) manifold
   int stable;		

   double h;	// small increment in the direction of v

   // Index specifying which iterate of fixed point we are interested in,
   // e.g. p_2 -> ifp==2
   int ifp;		
   	
   double p0[2];	// p0 = p+hv
   double p1[2];	// p1 = P(p0)

   // Linear segment approximating local invariant manifold
   double l[2*NPOINTS];
   double l4[DIM*NPOINTS];	// 4D version of l

   double err;		// error commited in approximating the manifold

   // Auxiliary variables
   int status, iter, i, j;
   double ti;
   char section_str[10];        // holds input string "SEC1", "SEC2" etc

   // 1. Input parameters from stdin.
   if(scanf("%le %s %le %d %le %le %le %le %le %d %d %le %d", 
	    &mu, section_str, &H, &k, p, p+1, v, v+1, &lambda, &n, &stable, &h, &ifp) < 13)
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

   // Estimate error commited in the linear approximation of the manifold
   err = err_mfld(mu,sec,H,k,p,v,lambda,stable,h);
   fprintf(stderr,"Estimated error of manifold: %le\n", err);

   // 2. Discretize the linear segment between $p0=p+hv$ and $p1=P(p0)$ into
   // a set of NPOINTS points.

   // We choose the small increment $h$ such that the estimated error commited
   // in the linear approximation of the manifold is smaller than 10^{-8}.
   //h=1.e-6;
   //fprintf(stderr, "h: %le\n", h);

   // Compute $p_0$
   p0[0] = p[0] + h*v[0]; 
   p0[1] = p[1] + h*v[1];

   // Compute $p_1$
   p1[0] = p0[0];
   p1[1] = p0[1];
   if(!stable) 	// unstable manifold
      status=prtbp_nl_2d(mu,sec,H,k,p1,&ti); 	// $p_1 = P(p_0)$
   else 	// stable manifold
      status=prtbp_nl_2d_inv(mu,sec,H,k,p1,&ti);	// $p_1 = P^{-1}(p_0)$
   if(status)
   {
      fprintf(stderr, "main: error computing Poincare map\n");
      return(1);
   }

   // Discretize linear segment
   disc(p0, p1, NPOINTS, l);

   // Lift points in linear segment from 2d to 4d
   status = lift(mu, sec, H, NPOINTS, l, l4);
   if(status)
   {
      fprintf(stderr, "main: error lifting point\n");
      return(1);
   }

   // 3. Iterate the (discretized) linear segment "n" times by the Poincare map,
   // i.e. compute its orbit (and print it to stdout). 
   for(iter=0;iter<n;iter++)
   {
	 if(!stable)	// unstable manifold
	 {
		for(i=0;i<NPOINTS;i++)
	    {
	       if(prtbp_nl(mu,sec,k,l4+DIM*i,&ti))
	       {
		      fprintf(stderr, "main: error computing Poincare map\n");
		      exit(EXIT_FAILURE);
	       }
	    }
	 }
	 else	// stable manifold
	 {
		for(i=0;i<NPOINTS;i++)
	    {
	       if(prtbp_nl_inv(mu,sec,k,l4+DIM*i,&ti))
	       {
		      fprintf(stderr, "main: error computing Poincare map\n");
		      exit(EXIT_FAILURE);
	       }
	    }
	 }

	 // Print iteration of linear segment
	 for(i=0;i<NPOINTS;i++)
	 {
		dblprint(l4+DIM*i, DIM);
		printf("\n");
	 }
	 printf("\n");
   }

   exit(EXIT_SUCCESS);
}
