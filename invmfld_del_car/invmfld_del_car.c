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
#include <rtbp.h>	// DIM
#include <hinv.h>	
#include <cardel.h>
#include <prtbp_nl_2d.h>
#include <prtbp_del_car.h>
#include <errmfld.h>
#include <disc.h>
#include <approxint_del_car.h>	// iterate_segment
#include <utils_module.h>		// dblcpy
#include "lift.h"

/**
  Main program.

  Input params (stdin): 
  
//    - mass parameter 
//    - Cartesian Poincare section "sec"
//    - energy value "H"
//    - number of cuts "k" with Poincare section
//    - fixed point "p"
//    - linear unstable (or stable) direction "v"
//    - unstable (or stable) eigenvalue "lambda"
//    - Delaunay Poincare section "sec_del"
//    - number of desired iterations by the Delaunay Poincare map "n"
//    - "stable": flag specifying whether to compute the unstable manifold
//    (unstable==0) or the stable manifold (stable==1).
//    - h: small increment in the direction of v

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
	/// Number of points in the discretization of unst segment.
	/// WARNING!!! Oddly, approxint_del_car does not work properly unless
	/// linear segment is discretized into very few points (e.g. 5).
	/// Probably, the more points we use, the higher the probability that
	/// prtbp_del_car fails.
	const int NPOINTS=11;

   double mu, H;
   section_t sec;
   int k;		// number of iterates of Poincare map

   double p[2];		// fixed point
   double v[2];		// stable/unstable direction
   double lambda;	// stable/unstable eigenvalue

   section_t sec_del;   // Delaunay Poincare section
   int n;		        // number of desired iterations of linear segment

   // "stable" flag specifies wheather we want to compute the unstable (=0)
   // or stable (=1) manifold
   int stable;		

   double h;	// small increment in the direction of v

   double p0[2];	// p0 = p+hv
   double p1[2];	// p1 = P(p0)

   // Linear segment approximating local invariant manifold
   double l[2*NPOINTS];

   // Linear segment approximating local invariant manifold (points in 4d)
   double l4_car[DIM*NPOINTS];

   // Linear segment approximating local invariant manifold (Delaunay coords)
   double l4_del[DIM*NPOINTS];

   double err;		// error commited in approximating the manifold

   // Auxiliary variables
   int status, iter, i, j;
   double ti;
   double l4_del_cpy[DIM*NPOINTS];  // copy of l4_del
   double l4_car_cpy[DIM*NPOINTS];  // copy of l4_car

   char section_str[10];        // holds input string "SEC1", "SEC2" etc
   char sec_del_str[10];        // holds input string "SECg", "SECg2" etc

   // 1. Input parameters from stdin.
   if(scanf("%le %s %le %d %le %le %le %le %le %s %d %d %le", 
        &mu, section_str, &H, &k, p, p+1, v, v+1, &lambda, sec_del_str, &n,
        &stable, &h) < 13)
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

   if (strcmp(sec_del_str,"SECg") == 0)
      sec_del = SECg;
   else if (strcmp(sec_del_str,"SECg2") == 0)
      sec_del = SECg2;
   else
   {
      perror("main: error reading sec_del string");
      exit(EXIT_FAILURE);
   }

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
   /*
   p1[0] = p0[0] + lambda*h*v[0];
   p1[1] = p0[1] + lambda*h*v[1];
   */

   // Discretize linear segment
   disc(p0, p1, NPOINTS, l);

   // IMPROVEMENT: instead of lift+cardel here, why not use cardel_2d???

   // Lift points in linear segment from \R^2 to \R^4
   status=lift(mu,sec,H,NPOINTS,l,l4_car);
   if(status)
   {
      perror("main: error lifting points in linear segment");
      exit(EXIT_FAILURE);
   }

   // Transform points in linear segment to Delaunay coordinates
   for(i=0; i<NPOINTS; i++)
       cardel(l4_car+DIM*i,l4_del+DIM*i);

   // Flow the (discretized) linear segment to Delaunay section before the
   // iteration.
   if(iterate_segment(mu,sec_del,1,1,l4_del,l4_car))
   {
	  fprintf(stderr,
			  "main: error flowing linear segment to Delaunay seciton"
			  " before the iteration\n");
      return(1);
   }

   dblcpy(l4_del_cpy, l4_del, DIM*NPOINTS);
   dblcpy(l4_car_cpy, l4_car, DIM*NPOINTS);

   // 3. Iterate the (discretized) linear segment "n" times 
   // by the Poincare map, i.e. compute its orbit (and print it to stdout). 
   for(iter=1;iter<=n;iter++)
   {
       dblcpy(l4_del, l4_del_cpy, DIM*NPOINTS);
       dblcpy(l4_car, l4_car_cpy, DIM*NPOINTS);

	    if(!stable)	// unstable manifold
	    {
		   if(iterate_segment(mu,sec_del,1,iter,l4_del,l4_car))
		   {
              fprintf(stderr, "main: error computing Poincare map\n");
              exit(EXIT_FAILURE);
		   }
	    }
	    else		// stable manifold
	    {
			/*
		   if(iterate_segment_inv(mu,sec_del,1,iter,l4_del,l4_car))
		   {
              fprintf(stderr, "main: error computing Poincare map\n");
              exit(EXIT_FAILURE);
		   }
		   */
	    }

		// Print iteration of linear segment
		for(i=0;i<NPOINTS;i++)
		{
		   if(printf("% .15le % .15le\n", 
					   l4_del[DIM*i+0], 
					   l4_del[DIM*i+1])<0)
		   {
			  perror("main: error writting output");
			  exit(EXIT_FAILURE);
		   }
		}
		//printf("\n");
   }

   // 4. Estimate error commited in the linear approximation of the manifold
   err = err_mfld(mu,sec,H,k,p,v,lambda,stable,h);
   fprintf(stderr,"Estimated error of manifold: %le\n", err);

   exit(EXIT_SUCCESS);
}
