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
#include <prtbp_2d.h>
#include <prtbp_del_car.h>
#include <errmfld.h>
#include <disc.h>	// disc

/// Number of points in discretization of linear segment
const int NPOINTS = 100; 

int lift(double mu, section_t sec, double H, int n, const double *l, 
        double *l4);

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
//    (unstable==0) or the stable manifold (stable==1).
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

   // Linear segment approximating local invariant manifold (points in 4d)
   double l4[DIM*NPOINTS];

   // Linear segment approximating local invariant manifold (Delaunay coords)
   double l4_del[DIM*NPOINTS];

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
      status=prtbp_2d(mu,sec,H,k,p1,&ti); 	// $p_1 = P(p_0)$
   else 	// stable manifold
      status=prtbp_2d_inv(mu,sec,H,k,p1,&ti);	// $p_1 = P^{-1}(p_0)$
   if(status)
   {
      fprintf(stderr, "main: error computing Poincare map\n");
      return(1);
   }

   // Discretize linear segment
   disc(p0, p1, NPOINTS, l);

   // Lift points in linear segment from \R^2 to \R^4
   status=lift(mu,sec,H,NPOINTS,l,l4);
   if(status)
   {
      perror("main: error lifting points in linear segment");
      exit(EXIT_FAILURE);
   }

   // Transform points in linear segment to Delaunay coordinates
   for(i=0; i<NPOINTS; i++)
       cardel(l4+DIM*i,l4_del+DIM*i);

   // 3. Iterate the (discretized) linear segment "n" times 
   // by the Poincare map, i.e. compute its orbit (and print it to stdout). 
   for(iter=0;iter<n;iter++)
   {
      for(j=1;j<=3;j++)
      {
	 for(i=0;i<NPOINTS;i++)
	 {
	    if(!stable)	// unstable manifold
	    {
	       if(prtbp_del_car(mu,SEC2,1,l4_del+DIM*i,l4+DIM*i,&ti))
	       {
              fprintf(stderr, "main: error computing Poincare map\n");
              exit(EXIT_FAILURE);
	       }
	    }
	    else		// stable manifold
	    {
	       if(prtbp_del_car_inv(mu,SEC2,1,l4_del+DIM*i,l4+DIM*i,&ti))
	       {
              fprintf(stderr, "main: error computing inverse Poincare map\n");
              exit(EXIT_FAILURE);
	       }
	    }
	 }

	 if(ifp==-1 || j==ifp)	// print only manifold of i-th fixed point p_i
	 {
	    // Print iteration of linear segment
	    for(i=0;i<NPOINTS;i++)
	    {
	       if(printf("% .15le % .15le\n", 
                       l4_del[DIM*i+2], 
                       l4_del[DIM*i+3])<0)
	       {
              perror("main: error writting output");
              exit(EXIT_FAILURE);
	       }
	    }
	    printf("\n");
	 }
      }
   }

   // 4. Estimate error commited in the linear approximation of the manifold
   err = err_mfld(mu,sec,H,k,p,v,lambda,stable,h);
   fprintf(stderr,"Estimated error of manifold: %le\n", err);

   exit(EXIT_SUCCESS);
}

int lift(double mu, section_t sec, double H, int n, const double *l, 
        double *l4)
{
    // Auxiliary variables
    int i, status;
    const double *p; 
    double *p4;

    for(i=0; i<n; i++)
    {
        p=l+2*i;
        p4=l4+DIM*i;

        p4[0]=p[0]; // x
        p4[1]=0;    // y
        p4[2]=p[1]; // p_x
        status=hinv(mu,sec,H,p4);
       if(status)
       {
          fprintf(stderr, "lift: error lifting point\n");
          return(1);
       }
    }
    return(0);
}
