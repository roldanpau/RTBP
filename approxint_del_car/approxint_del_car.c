/*! \file
    \brief Approximate Intersection of Invariant Manifolds

    $Author: roldan $
    $Date: 2013-03-26 22:10:02 $
*/

#include <stdio.h>
#include <stdlib.h>     // EXIT_SUCCESS, EXIT_FAILURE
#include <stdbool.h>	// bool data type
#include <string.h>	// memcpy
#include <math.h>	// atan2
#include <assert.h>
#include <hinv.h>
#include <cardel.h>
#include <prtbp_2d.h>
#include <prtbp_del_car.h>
#include <disc.h>	// disc
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
const int NPOINTS = 5;

int 
u_i (double mu, int k, double a, double *l4_del, double *l4, int *idx);
int 
s_i (double mu, double H, int k, double z[2], double a, double *l, int *idx);

// Obs! parameter lambda_u is not used?

int 
approxint_del_car_unst (double mu, double H, int k, double p[2], double v[2],
      double lambda, double h, double a, 
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
   int status, iter, i;
   double ti;
   double l_bak[2*NPOINTS];    // Aux copy of l

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
   status=prtbp_2d(mu,SEC2,H,k,p1,&ti);        // $p_1 = P(p_0)$
   if(status)
   {
      fprintf(stderr, "approxint_del_car: error computing Poincare map\n");
      return(1);
   }

   // Discretize linear segment
   disc(p0, p1, NPOINTS, l);

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

   // Make a copy of l in order to not modify it.
   memcpy(l_bak, l, 2*NPOINTS*sizeof(double));

   // Iterate the (discretized) linear segment twice until we are 
   // at the right "eye" of the resonance.
   for(i=0;i<NPOINTS;i++)
   {
       if(prtbp_del_car(mu,SEC2,3,l4_del+DIM*i,l4+DIM*i,&ti))
       {
          fprintf(stderr, "approxint_del_car: error computing Poincare map\n");
          return(1);
       }
   }

   // Iterate the (discretized) linear segment "iter" times by the Poincare map
   for(iter=1;iter<=MAXITER;iter++)
   {
      status=u_i(mu, 3, a, l4_del, l4, &i);
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
      // Manifold does not intersect $G$ axis!!
      return(2);
   }

   // Manifold intersects $G$ axis.
   *piter = iter;

   // Endpoints of segment u_i: 
   // P[2] = (x, p_x) = (l_bak[2i], l_bak[2i+1]),
   // Q[2] = (x, p_x) = (l_bak[2i+2], l_bak[2i+3]).
   *h_1 = (l_bak[2*i+1]-p[1])/v[1];
   *h_2 = (l_bak[2*i+3]-p[1])/v[1];

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
   status=prtbp_2d_inv(mu,SEC2,H,k,p1,&ti);        // $p_1 = \sixmap^{-1}(p_0)$
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
      status=s_i(mu, H, k, z, a, l, &i);
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
// H
//    energy value
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
u_i (double mu, int k, double a, double *l4_del, double *l4, int *idx)
{
   // Auxiliary variables
   int status, iter, i;
   double ti;
   double dx,dy;

   // Approximate splitting half-angle
   //double alpha2;

   // 3. Iterate the (discretized) linear segment one more time by the
   // Poincare map
   for(i=0;i<NPOINTS;i++)
   {
       if(prtbp_del_car(mu,SEC2,k,l4_del+DIM*i,l4+DIM*i,&ti))
       {
          fprintf(stderr, "u_i: error computing Poincare map of %d-th point\n", i);
          return(1);
       }
   }

   // We look for the first unst segment U_i that crosses the $G$ axis.
   for(i=0; i<(NPOINTS-1); i++)
   {
      // Endpoints of segment U_i: 
      // P[2] = (g, G) = (l4_del[4i+2], l4_del[4i+3]),
      // Q[2] = (g, G) = (l4_del[4(i+1)+2], l4_del[4(i+1)+3]).
  
       // Since dg/dt<0, it is enough to check if the angle
       // has been reset from 0 to 2\pi
      if(l4_del[DIM*(i+1)+2] < l4_del[DIM*i+2]) break;
   }

   if(i==(NPOINTS-1))
   {
      // Manifold does not intersect $G$ axis
      *idx=-1;
      return(0);
   }

   // Manifold intersects the $G$ axis.

   //dx = l[2*i+1]-l[2*i+3];	// P2 - Q2
   //dy = l[2*i]-l[2*i+2];	// P1 - Q1
   //alpha2 = atan2(dy,dx);	
   //fprintf(stderr, "Approximate splitting half-angle: %f\n", alpha2);

   fprintf(stderr, "Approximate interval: %d\n", i);
   *idx=i;
   return(0);
}

/** 
  \note We have problems computing some invariant manifolds, because suddenly
  they look "broken". To prevent this, we check that the manifolds are indeed
  continuous, i.e. we check that each two consecutive points in the manifold
  are close together (closer than say 0.1).
  */ 

int 
s_i (double mu, double H, int k, double z[2], double a, double *l, int *idx)
{
   // Auxiliary variables
   int status, iter, i;
   double ti;
   double dx,dy;

   // Approximate splitting half-angle
   double alpha2;

   // 3. Iterate the (discretized) linear segment one more time by the
   // Poincare map
   for(i=0;i<NPOINTS;i++)
   {
	 if(prtbp_2d_inv(mu,SEC2,H,k,l+2*i,&ti))
	 {
	    fprintf(stderr, "s_i: error computing Poincare map\n");
	    return(1);
	 }
   }

   // We look for the first st segment S_i that crosses the $x$ axis.
   for(i=0; i<(NPOINTS-1); i++)
   {
      // Endpoints of segment S_i: 
      // P[2] = (x, p_x) = (l[2i], l[2i+1]),
      // Q[2] = (x, p_x) = (l[2i+2], l[2i+3]).

      // Check that the manifolds are indeed continuous, i.e. we check that
      // each two consecutive points in the manifold  are close together.
      if( fabs(l[2*i+2]-l[2*i]) + fabs(l[2*i+3]-l[2*i+1]) > 0.1 )
      {
	 // Segment S_i is too large
	 fprintf(stderr, "s_i: segment is too large!\n");
	 return(1);
      }

      if(((l[2*i+1]-a)*(l[2*i+3]-a)<=0) && 
	    fabs(l[2*i]-z[0]) < 0.02)	// point belongs to primary family
	 break;
   }

   if(i==(NPOINTS-1))
   {
      // Manifold does not intersect $x$ axis
      *idx=-1;
      return(0);
   }

   // Manifold intersects the $x$ axis.

   //dx = l[2*i+3]-l[2*i+1];	// Q2 - P2
   //dy = l[2*i+2]-l[2*i];	// Q1 - P1
   //alpha2 = atan2(dy,dx);	
   //fprintf(stderr, "Approximate splitting half-angle: %f\n", alpha2);

   fprintf(stderr, "Approximate interval: %d\n", i);
   *idx=i;
   return(0);
}
