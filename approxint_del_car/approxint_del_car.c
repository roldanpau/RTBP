/*! \file
    \brief Approximate Intersection of Invariant Manifolds

    $Author: roldan $
    $Date: 2013-03-26 22:10:02 $
*/

#include <stdio.h>
#include <stdlib.h>     // EXIT_SUCCESS, EXIT_FAILURE
#include <stdbool.h>	// bool data type
#include <string.h>	    // memcpy
#include <math.h>	    // atan2
#include <assert.h>
#include <hinv.h>
#include <cardel.h>
#include <section.h>    // section_t
#include <prtbp_nl_2d.h>
#include <prtbp_del_car.h>
#include <disc.h>	   
#include <lift.h>

/*! \brief 
  Max number of iterations of unst segment before we give up looking for
  intersection.
 */
const int MAXITER = 100; 

/// Number of points in the discretization of unst segment.
/// WARNING!!! Oddly, approxint_del_car does not work properly unless 
/// linear segment is discretized into very few points (e.g. 5).
/// Probably, the more points we use, the higher the probability that 
/// prtbp_del_car fails.
const int NPOINTS = 50;

int 
u_i (double mu, section_t sec, int k, double a, double *l4_del, double *l4, 
        int *idx);
int 
s_i (double mu, double H, int k, double z[2], double a, double *l, int *idx);

// Obs! parameter lambda_u is not used?

int 
approxint_del_car_unst (double mu, section_t sec, double H, int k, 
        double p[2], double v[2], double lambda, double h, double a, 
        int *piter, double *h_1, double *h_2, double z[2])
{
   double p0[2];        // p0 = p+hv
   double p1[2];        // p1 = P(p0)

   // Linear segment approximating local invariant manifold
   double l[2*NPOINTS];

   // Linear segment approximating local invariant manifold (points in 4d)
   double l4[DIM*NPOINTS];

   // Linear segment approximating local invariant manifold (Delaunay coords)
   double l4_del[DIM*NPOINTS];

   // Auxiliary variables
   int status, iter, i, iskip;
   double ti;

   // we DO NOT assume that $v=(x,p_x)$ points "to the right", i.e.
   // we DO NOT assume that the first component of $v$ is $x>0$
   //assert(v[0]>0);

   // 2. Discretize the linear segment between $p0=p+hv$ and $p1=P(p0)$ into
   // a set of NPOINTS points.

   // Compute $p_0$
   p0[0] = p[0] + h*v[0];
   p0[1] = p[1] + h*v[1];

   // Compute $p_1$
   p1[0] = p0[0];
   p1[1] = p0[1];
   status=prtbp_nl_2d(mu,SEC2,H,k,p1,&ti);        // $p_1 = P(p_0)$
   if(status)
   {
      fprintf(stderr, "approxint_del_car: error computing Poincare map\n");
      return(1);
   }

   // Discretize linear segment
   disc(p0, p1, NPOINTS, l);

   // IMPROVEMENT: Use cardel_2d instead of lift+cardel.

   // Lift points in linear segment from \R^2 to \R^4
   status=lift(mu,SEC2,H,NPOINTS,l,l4);
   if(status)
   {
      perror("main: error lifting points in linear segment");
      exit(EXIT_FAILURE);
   }

   // Transform points in linear segment to Delaunay coordinates
   for(i=0; i<NPOINTS; i++)
       cardel(l4+DIM*i,l4_del+DIM*i);

   // Flow the (discretized) linear segment to Delaunay section before the
   // iteration.
   status=u_i(mu, sec, 1, a, l4_del, l4, &i);

   // Iterate the (discretized) linear segment "iter" times by the Poincare map
   for(iter=1;iter<=MAXITER;iter++)
   {
      //fprintf(stderr, "DEBUG: before %d iteration of linear segment: l=%.15le, g=%.15le\n", iter, l4_del[0], l4_del[2]);
      status=u_i(mu, sec, 1, a, l4_del, l4, &i);
      //fprintf(stderr, "DEBUG: after %d iteration of linear segment: l=%.15le, g=%.15le\n", iter, l4_del[0], l4_del[2]);
      if(status)
      {
          fprintf(stderr, 
                  "approxint_del_car: error during %d-th iteration of linear segment\n", iter);
          return(1);
      }
      if(i>=0)	// intersection found
	 break;
   }

   if(iter==(MAXITER+1))
   {
      // Manifold does not intersect line $g=a$!!
      return(2);
   }

   // Manifold intersects line $g=a$
   *piter = iter;

   // Endpoints of segment u_i: 
   // P[2] = (x, p_x) = (l[2i], l[2i+1]),
   // Q[2] = (x, p_x) = (l[2i+2], l[2i+3]).
   *h_1 = (l[2*i+1]-p[1])/v[1];
   *h_2 = (l[2*i+3]-p[1])/v[1];

   // return approximate intersection point
   z[0] = l4_del[DIM*i+2];
   z[1] = l4_del[DIM*i+3];
   return(0);
}

int 
approxint_del_car_st (double mu, double H, int k, double p[2], double v[2],
      double lambda, double h, double a, 
      int *piter, double *h_1, double *h_2, double z[2])
{
   double p0[2];        // p0 = p+hv
   double p1[2];        // p1 = P(p0)

   // Linear segment approximating local invariant manifold
   double l[2*NPOINTS];

   // Auxiliary variables
   int status, iter, i;
   double ti;
   double l_bak[2*NPOINTS];    // Aux copy of l

   fprintf(stderr, "approxint_del_car_st: this function is not implemented!\n");
   exit(1);

   // we DO NOT assume that $v=(x,p_x)$ points "to the right", i.e.
   // we DO NOT assume that the first component of $v$ is $x>0$
   //assert(v[0]>0);

   // 2. Discretize the linear segment between $p0=p+hv$ and $p1=P(p0)$ into
   // a set of NPOINTS points.

   // Compute $p_0$
   p0[0] = p[0] + h*v[0];
   p0[1] = p[1] + h*v[1];

   // Compute $p_1$
   p1[0] = p0[0];
   p1[1] = p0[1];
   status=prtbp_nl_2d_inv(mu,SEC2,H,k,p1,&ti);        // $p_1 = \sixmap^{-1}(p_0)$
   if(status)
   {
      fprintf(stderr, "approxint_del_car: error computing Poincare map\n");
      return(1);
   }

   // Discretize linear segment
   disc(p0, p1, NPOINTS, l);

   // Make a copy of l in order to not modify it.
   memcpy(l_bak, l, 2*NPOINTS*sizeof(double));

   // Iterate the (discretized) linear segment "iter" times by the Poincare map
   for(iter=1;iter<=MAXITER;iter++)
   {
      //status=s_i(mu, H, k, z, a, l, &i);
      if(status)
      {
	 fprintf(stderr, 
	       "approxint_del_car: error during %d-th iteration of linear segment\n",
	       iter);
	 return(1);
      }
      if(i>=0)	// intersection found
	 break;
   }

   if(iter==(MAXITER+1))
   {
      // Manifold does not intersect $x$ axis!!
      return(2);
   }

   // Manifold intersects $x$ axis.
   *piter = iter;

   // Endpoints of segment s_i: 
   // P[2] = (x, p_x) = (l_bak[2i], l_bak[2i+1]),
   // Q[2] = (x, p_x) = (l_bak[2i+2], l_bak[2i+3]).
   *h_1 = (l_bak[2*i+1]-p[1])/v[1];
   *h_2 = (l_bak[2*i+3]-p[1])/v[1];

   // update approximate intersection point
   z[0] = l[2*i];
   z[1] = l[2*i+1];
   return(0);
}

// name OF FUNCTION: u_i
//
// PURPOSE
// =======
// This function checks whether the n-th iteration does intersect the $x$
// axis or not. 
// If it does intersect the $x$ axis, it returns the unstable segment u_i
// such that its image $U_i$ crosses the $x$ axis.
//
// We discretize the unst domain into a set of NPOINTS segments u_1, u_2,
// ..., u_NPOINTS. 
// The image under iteration of this discretized version of the unst manifold
// is a set of segments U_1, U_2, ..., U_NPOINTS that approximates the
// nonlinear unst manifold.
// 
// We look for the first unst segment U_i that intersects the $x$ axis.
// Therefore, U_i contains an intersection point, and the segment u_i in the
// fundamental domain contains the preimage of an intersection point.
// It returns the unstable segment u_i such that its image $U_i$ crosses the
// $x$ axis.
//
// PARAMETERS
// ==========
// mu
//    mass parameter for the RTBP
// sec        
//    Poincare section: sec={SEC1,SEC2}
// k
//    number of cuts with poincare section per iterated poincare map
//    $\sixmap$.
// a
//    axis $p_x=a$ parallel to the $x$ axis.
// lu
//    NPOINTS points approximating linear unstable fundamental domain.
//    On exit, this vector is modified with a new iteration of the
//    fundamental domain.
// idx
//    On exit, it contains the index $i$ corresponding to interval $u_i$
//    such that its image $U_i$ crosses the $x$ axis.
//    If the manifold does not cross the axis, we return i=-1.
// 
// RETURN VALUE
// ============
// Returns a non-zero error code to indicate an error and 0 to indicate
// success:
//    1: Problems computing the Poincare iterates.

int 
u_i (double mu, section_t sec, int k, double a, double *l4_del, double *l4, 
        int *idx)
{
   // Auxiliary variables
   int status, iter, i;
   double ti;
   double dx,dy;
   double g1, g2;
   double l1, l2;

   // Approximate splitting half-angle
   //double alpha2;

   // 3. Iterate the (discretized) linear segment one more time by the
   // Poincare map
   for(i=0;i<NPOINTS;i++)
   {
       if(prtbp_del_car(mu,sec,k,l4_del+DIM*i,l4+DIM*i,&ti))
       {
          fprintf(stderr, "u_i: error computing Poincare map of %d-th point\n", i);
          return(1);
       }
   }

   // We look for the first unst segment U_i that crosses the line $g=a$.
   for(i=0; i<(NPOINTS-1); i++)
   {
       l1=l4_del[DIM*i];
       l2=l4_del[DIM*(i+1)];
  
      // Endpoints of segment U_i: 
      // P[2] = (g_1, G_1) = (l4_del[4i+2], l4_del[4i+3]),
      // Q[2] = (g_2, G_2) = (l4_del[4(i+1)+2], l4_del[4(i+1)+3]).
       g1=l4_del[DIM*i+2];
       g2=l4_del[DIM*(i+1)+2];

       if(sec==SEC1)
       {
          // For the upper branch, since dg/dt<0, it is enough to check 
          // if the angle has passed from >a to <a
          if(g1 > a && g2 < a) break;
          
          // For the lower branch, since dg/dt>0, it is enough to check 
          // if the angle has passed from <a to >a
          //if(g1 < a && g2 > a) break;
       }
       else if(sec==SEC2)
       {
          // For the upper branch, since dg/dt<0, it is enough to check 
          // if the angle has been reset from 0 to 2\pi
          if(g2 > g1) {
              //fprintf(stderr, "DEBUG: g=%.15le, gprime=%.15le\n", l4_del[DIM*i+2], l4_del[DIM*(i+1)+2]);
              break;
          }
          
          // For the lower branch, since dg/dt>0, it is enough to check 
          // if the angle has been reset from 2\pi to 0
          //if(g2 < g1) break;
       }
       else if(sec==SECg)
       {
          // For the upper branch, since dl/dt<0, it is enough to check 
          // if the angle has been reset from -\pi to \pi
          if(l2 > l1) {
              //fprintf(stderr, "DEBUG: g=%.15le, gprime=%.15le\n", l4_del[DIM*i+2], l4_del[DIM*(i+1)+2]);
              break;
          }
          
          // For the lower branch, since dl/dt>0, it is enough to check 
          // if the angle has been reset from \pi to -\pi
          //if(l2 < l1) break;
       }
       else
       {
          fprintf(stderr, "u_i: unknown section type. Exiting\n");
          exit(1);
       }
   }

   if(i==(NPOINTS-1))
   {
      // Manifold does not intersect line $g=a$
      *idx=-1;
      return(0);
   }

   // Manifold intersects the line $g=a$.

   //dx = l[2*i+1]-l[2*i+3];	// P2 - Q2
   //dy = l[2*i]-l[2*i+2];	// P1 - Q1
   //alpha2 = atan2(dy,dx);	
   //fprintf(stderr, "Approximate splitting half-angle: %f\n", alpha2);

   fprintf(stderr, "Approximate interval: %d\n", i);
   *idx=i;
   return(0);
}
