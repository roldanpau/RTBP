// ==========================================================
// Invariant Manifolds of fixed point of Poincare map in RTBP
// ==========================================================
// FILE:          $RCSfile: invmflddel.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2012-06-06 12:40:08 $
//
// PURPOSE
// =======
// Consider the Restricted Three Body Problem. Assume that energy is fixed:
// $H(l,L,g,G)=\bar H$.
// Let "sec" be a Poincare section, and $P(g,G)$ be the associated 2D
// Poincare map (consisting of $k$ cuts with the Poincare section).
// Let $p(g,G)$ be a fixed point of the Poincare map associated to a
// periodic orbit of the flow.
// This function computes the unstable/stable manifolds $W^u(p), W^s(p)$ of
// the point $p$. A manifold is represented numerically as a sets of points.
//
// OVERALL METHOD
// ==============
//
// 1. Input parameters from stdin:
// 
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
//
// 2. Let $h$ be a small increment. Discretize the linear segment between
// $p+hv_u$ and $P(p+hv_u)$ into a set of NPOINTS points.
//
// 3. Iterate the (discretized) linear segment "n" times by the Poincare map,
// i.e. compute its orbit (and print it to stdout). 
//
// 4. Estimate error commited in the linear approximation of the manifold

// NOTES
// =====
// If the flag "stable" specifies the unstable manifold (0), we iterate the
// forward Poincare map $P$.
// If the flag "stable" specifies the stable manifold (1), we iterate the
// backward Poincare map $P^{-1}$.
//
// OBSOLETE:
// We take the small increment $h$ as a function of the stable/unstable
// eigenvalue. The strongest the contraction/expansion is, the smallest h we
// take.
// In particular, in the case of the stable direction, we take 
//    h=1.e-5*lambda,
// where \lambda is the stable eigenvalue. In the case of the unst direction,
//    h=1.e-5/lambda,
// where \lambda is the unstable eigenvalue.
//
// NEW:
// We choose the small increment $h$ such that the estimated error commited
// in the linear approximation of the manifold is smaller than 10^{-8}.

#include <stdio.h>
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE
#include <math.h>	// sqrt
#include <prtbpdel.h>	// section_t
#include <prtbpdel_2d.h>	// prtbp_del_2d, prtbp_del_2d_inv
#include <disc.h>	// disc

// Number of points in discretization of linear segment
const int NPOINTS = 100; 

int err_mfld(double mu, section_t sec, double H, int k, double p[2], 
      double v[2], double lambda, int stable, double h);

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
      status=prtbp_del_2d(mu,sec,H,k,p1,&ti); 	// $p_1 = P(p_0)$
   else 	// stable manifold
      status=prtbp_del_2d_inv(mu,sec,H,k,p1,&ti);	// $p_1 = P^{-1}(p_0)$
   if(status)
   {
      fprintf(stderr, "main: error computing Poincare map\n");
      return(1);
   }

   // Discretize linear segment
   disc(p0, p1, NPOINTS, l);

   // 3. Iterate the (discretized) linear segment "n" times by the Poincare map,
   // i.e. compute its orbit (and print it to stdout). 
   for(iter=0;iter<n;iter++)
   {
      for(j=1;j<=k;j++)
      {
	 for(i=0;i<NPOINTS;i++)
	 {
	    if(!stable)	// unstable manifold
	    {
	       if(prtbp_del_2d(mu,sec,H,1,l+2*i,&ti))
	       {
		  fprintf(stderr, "main: error computing Poincare map\n");
		  exit(EXIT_FAILURE);
	       }
	    }
	    else		// stable manifold
	    {
	       if(prtbp_del_2d_inv(mu,sec,H,1,l+2*i,&ti))
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
	       if(printf("% .15le % .15le\n", l[2*i], l[2*i+1])<0)
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
   err_mfld(mu,sec,H,k,p,v,lambda,stable,h);

   exit(EXIT_SUCCESS);
}

int err_mfld(double mu, section_t sec, double H, int k, double p[2], double v[2], 
      double lambda, int stable, double h)
{
   double p1[2];	// p1 = P(p+hv)
   double e[2];
   double err;

   // Auxiliary variables
   int status;
   double ti;

   // Compute $p_1$
   p1[0] = p[0] + h*v[0];
   p1[1] = p[1] + h*v[1];
   if(!stable) 	// unstable manifold
      status=prtbp_del_2d(mu,sec,H,k,p1,&ti); 	// $p_1 = P(p_0)$
   else 	// stable manifold
      status=prtbp_del_2d_inv(mu,sec,H,k,p1,&ti);	// $p_1 = P^{-1}(p_0)$
   if(status)
   {
      fprintf(stderr, "main: error computing Poincare map\n");
      return(1);
   }

   if(!stable) // unstable manifold
   {
      e[0] = p1[0]-p[0]-lambda*h*v[0];
      e[1] = p1[1]-p[1]-lambda*h*v[1];
   }
   else
   {
      e[0] = p1[0]-p[0]-1.0/lambda*h*v[0];
      e[1] = p1[1]-p[1]-1.0/lambda*h*v[1];
   }
   err = sqrt(e[0]*e[0]+e[1]*e[1]);
   fprintf(stderr,"Estimated error of manifold: %le\n", err);
}
