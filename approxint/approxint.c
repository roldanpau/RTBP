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
#include <rtbp.h>	// DIM
#include <lift.h>
#include <assert.h>
#include <prtbp_nl_2d_module.h>	// prtbp_nl_2d, prtbp_nl_2d_inv
#include <prtbp_nl.h>	// prtbp_nl, prtbp_nl_inv
#include <disc.h>	// disc
#include <utils_module.h>	// l2_norm

/*! \brief 
  Max number of iterations of unst segment before we give up looking for
  intersection.
 */
const int MAXITER = 100; 

/// Number of points in the discretization of unst segment.
const int NPOINTS = 100;

int 
u_i (double mu, double H, int k, double z[2], double a, double *l, int *idx);
int 
s_i (double mu, double H, int k, double z[2], double a, double *l, int *idx);

// NOTES
// =====
// Mental note: An alternative way to check that intersection belongs to
// primary family may be to check that piter is always increasing. If it was
// not, this would mean that the manifolds have developed a loop and a
// secondary intersection has appeared.
//
// Obs! parameter lambda_u is not used?

int 
approxint_unst (double mu, double H, int k, double p[2], double v[2],
      double lambda, double h, double a, 
      int *piter, double *h_1, double *h_2, double z[2])
{
   double p0[2];        // p0 = p+hv
   double p1[2];        // p1 = P(p0)

   // Linear segment approximating local invariant manifold
   double l[2*NPOINTS];
   double l4[DIM*NPOINTS];  // 4D version of l

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
   status=prtbp_nl_2d(mu,SEC2,H,k,p1,&ti);        // $p_1 = P(p_0)$
   if(status)
   {
      fprintf(stderr, "approxint: error computing Poincare map\n");
      return(1);
   }

   // Discretize linear segment
   disc(p0, p1, NPOINTS, l);

   // Lift points in linear segment from 2d to 4d
   status = lift(mu, SEC2, H, NPOINTS, l, l4);
   if(status)
   {
      fprintf(stderr, "main: error lifting point\n");
      return(1);
   }

   // Make a copy of l in order to not modify it.
   memcpy(l_bak, l, 2*NPOINTS*sizeof(double));

   // Iterate the (discretized) linear segment "iter" times by the Poincare map
   for(iter=1;iter<=MAXITER;iter++)
   {
      status=u_i(mu, H, k, z, a, l4, &i);
      if(status)
      {
	 fprintf(stderr, 
	       "approxint: error during %d-th iteration of linear segment\n",
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

   // Endpoints of segment u_i: 
   // P[2] = (x, p_x) = (l_bak[2i], l_bak[2i+1]),
   // Q[2] = (x, p_x) = (l_bak[2i+2], l_bak[2i+3]).
   *h_1 = (l_bak[2*i+1]-p[1])/v[1];
   *h_2 = (l_bak[2*i+3]-p[1])/v[1];

   // update approximate intersection point
   z[0] = l4[DIM*i];
   z[1] = l4[DIM*i+2];
   return(0);
}

// NOTES
// =====
// Mental note: An alternative way to check that intersection belongs to
// primary family may be to check that piter is always increasing. If it was
// not, this would mean that the manifolds have developed a loop and a
// secondary intersection has appeared.
//
// Obs! parameter lambda_u is not used?

int 
approxint_st (double mu, double H, int k, double p[2], double v[2],
      double lambda, double h, double a, 
      int *piter, double *h_1, double *h_2, double z[2])
{
   double p0[2];        // p0 = p+hv
   double p1[2];        // p1 = P(p0)

   // Linear segment approximating local invariant manifold
   double l[2*NPOINTS];
   double l4[DIM*NPOINTS];  // 4D version of l

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
   status=prtbp_nl_2d_inv(mu,SEC2,H,k,p1,&ti);        // $p_1 = P(p_0)$
   if(status)
   {
      fprintf(stderr, "approxint: error computing Poincare map\n");
      return(1);
   }

   // Discretize linear segment
   disc(p0, p1, NPOINTS, l);

   // Lift points in linear segment from 2d to 4d
   status = lift(mu, SEC2, H, NPOINTS, l, l4);
   if(status)
   {
      fprintf(stderr, "main: error lifting point\n");
      return(1);
   }

   // Make a copy of l in order to not modify it.
   memcpy(l_bak, l, 2*NPOINTS*sizeof(double));

   // Iterate the (discretized) linear segment "iter" times by the Poincare map
   for(iter=1;iter<=MAXITER;iter++)
   {
      status=s_i(mu, H, k, z, a, l4, &i);
      if(status)
      {
	 fprintf(stderr, 
	       "approxint: error during %d-th iteration of linear segment\n",
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

   // Endpoints of segment u_i: 
   // P[2] = (x, p_x) = (l_bak[2i], l_bak[2i+1]),
   // Q[2] = (x, p_x) = (l_bak[2i+2], l_bak[2i+3]).
   *h_1 = (l_bak[2*i+1]-p[1])/v[1];
   *h_2 = (l_bak[2*i+3]-p[1])/v[1];

   // update approximate intersection point
   z[0] = l4[DIM*i];
   z[1] = l4[DIM*i+2];
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
// z[2]
//    (input) contains the previous intersection point (for the previous
//    energy level).
// a
//    axis $p_x=a$ parallel to the $x$ axis.
// l4
//    NPOINTS 4d points approximating linear unstable fundamental domain.
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
u_i (double mu, double H, int k, double z[2], double a, double *l4, int *idx)
{
   // Auxiliary variables
   int status, iter, i;
   double ti;
   double dx,dy;

   // difference between two consecutive points on the manifold
   double v[2];	

   // Approximate splitting half-angle
   //double alpha2;

   // 3. Iterate the (discretized) linear segment one more time by the
   // Poincare map
   for(i=0;i<NPOINTS;i++)
   {
	 if(prtbp_nl(mu,SEC2,k,l4+DIM*i,&ti))
	 {
	    fprintf(stderr, "u_i: error computing Poincare map\n");
	    return(1);
	 }
   }

   // We look for the first unst segment U_i that crosses the $p_x=0$ axis.
   for(i=0; i<(NPOINTS-1); i++)
   {
      // Endpoints of segment U_i: 
      // P[2] = (x, p_x) = (l4[4i], l4[4i+2]),
      // Q[2] = (x, p_x) = (l4[4(i+1)], l4[4(i+1)+2]).

	   v[0] = l4[DIM*i]-l4[DIM*(i+1)];
	   v[1] = l4[DIM*i+2]-l4[DIM*(i+1)+2];
  
      if(((l4[DIM*i+2]-a)*(l4[DIM*(i+1)+2]-a)<=0) && 
	    l2_norm(v, 2) < 0.5)	// consecutive points are "close enough"
	 break;
   }

   if(i==(NPOINTS-1))
   {
      // Manifold does not intersect $x$ axis
      *idx=-1;
      return(0);
   }

   // Manifold intersects the $x$ axis.

   //dx = l[2*i+1]-l[2*i+3];	// P2 - Q2
   //dy = l[2*i]-l[2*i+2];	// P1 - Q1
   //alpha2 = atan2(dy,dx);	
   //fprintf(stderr, "Approximate splitting half-angle: %f\n", alpha2);

   fprintf(stderr, "Approximate interval: %d\n", i);
   *idx=i;
   return(0);
}

int 
s_i (double mu, double H, int k, double z[2], double a, double *l4, int *idx)
{
   // Auxiliary variables
   int status, iter, i;
   double ti;
   double dx,dy;

   // difference between two consecutive points on the manifold
   double v[2];	

   // Approximate splitting half-angle
   //double alpha2;

   // 3. Iterate the (discretized) linear segment one more time by the
   // Poincare map
   for(i=0;i<NPOINTS;i++)
   {
	 if(prtbp_nl_inv(mu,SEC2,k,l4+DIM*i,&ti))
	 {
	    fprintf(stderr, "u_i: error computing Poincare map\n");
	    return(1);
	 }
   }

   // We look for the first st segment S_i that crosses the $p_x=0$ axis.
   for(i=0; i<(NPOINTS-1); i++)
   {
      // Endpoints of segment S_i: 
      // P[2] = (x, p_x) = (l4[4i], l4[4i+2]),
      // Q[2] = (x, p_x) = (l4[4(i+1)], l4[4(i+1)+2]).

	   v[0] = l4[DIM*i]-l4[DIM*(i+1)];
	   v[1] = l4[DIM*i+2]-l4[DIM*(i+1)+2];
  
      if(((l4[DIM*i+2]-a)*(l4[DIM*(i+1)+2]-a)<=0) && 
	    l2_norm(v, 2) < 0.5)	// consecutive points are "close enough"
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

