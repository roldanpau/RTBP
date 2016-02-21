// ==========================================
// Approx Intersection of Invariant Manifolds
// ==========================================
// FILE:          $RCSfile: approxintdel.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2012-06-06 10:15:20 $
//
// FUNCTIONS
// =========
//
// approxint_del_unst
// approxint_del_st
// ------------------
// Consider the Restricted Three Body Problem. Assume that energy is fixed:
// $H(l,L,g,G)=\bar H$.
// Let "sec" be a Poincare section, and $P(g,G)$ be the associated 2D
// Poincare map.
// Let $p(g,G)$ be a fixed point of the Poincare map associated to the 3:1
// resonant periodic orbit of the flow.
// Let $W^u(p)$ be the unstable manifold of $p$.
// Let $\{g=a\}$ be a line parallel to the $G$ axis.
// This function computes an approximation to the "first" intersection of the
// unstable manifold with the line $g=a$.
// It returns the unstable segment u_i with endpoints $(h_1, h_2)$ containing
// the approximate root $p_u=p + h_u v_u$ with $h_u \in (h_1,h_2)$.
//
// u_i_del
// s_i_del
// -------
// This function checks whether the n-th iteration does intersect the $G$
// axis or not. 
// If it does intersect the $G$ axis, it returns the unstable segment u_i
// such that its image $U_i$ crosses the $G$ axis.

#include <stdio.h>
#include <stdlib.h>     // EXIT_SUCCESS, EXIT_FAILURE
#include <stdbool.h>	// bool data type
#include <string.h>	// memcpy
#include <math.h>	// atan2

#include <prtbpdel.h>		// SEC2
#include <prtbpdel_2d.h>	// prtbp_del_2d, prtbp_del_2d_inv
#include <disc.h>		// disc

// Max number of iterations of unst segment before we give up looking for
// intersection.
const int MAXITER = 40; 

// Number of points in the discretization of unst segment.
const int NPOINTS = 100;

int 
u_i_del (double mu, double H, int k, double z[2], double a, double *l, int *idx);
int 
s_i_del (double mu, double H, int k, double z[2], double a, double *l, int *idx);

// name OF FUNCTION: approxint_del_unst
//
// PURPOSE
// =======
// Consider the Restricted Three Body Problem. Assume that energy is fixed:
// $H(l,L,g,G)=\bar H$.
// Let SEC2: {l=pi} be a Poincare section, and $P(g,G)$ be the associated 2D
// Poincare map.
// Let $p(g,G)$ be a fixed point of the Poincare map associated to the 3:1
// resonant periodic orbit of the flow.
// Let $W^u(p)$ be the unstable manifold of $p$.
// Let $\{g=a\}$ be a line parallel to the $G$ axis.
// This function computes an approximation to the "first" intersection of the
// unstable manifold with the line $g=a$.
// It returns the unstable segment u_i in the form $(h_1, h_2)$ containing
// the approximate root $p_u=p + h_u v_u$ with $h_u \in (h_1,h_2)$.
// 
// We consider the unstable fundamental segment between the two points 
// $p+h_u v_u$ and $P(p+h_u v_u)$.
// We iterate the unstable fundamental domain until it intersects the $G$
// axis.
//
// We discretize the unst domain into a set of NPOINTS segments u_1, u_2,
// ..., u_NPOINTS. 
// The image under iteration of this discretized version of the unst manifold
// is a set of segments U_1, U_2, ..., U_NPOINTS that approximates the
// nonlinear unst manifold.
// 
// We look for the first unst segment U_i that intersects the $G$ axis.
// Therefore, U_i contains an intersection point, and the segment u_i in the
// fundamental domains contain the preimage of an intersection point.
// It returns the unstable segment u_i with endpoints $(h_1, h_2)$ containing
// the approximate root $p_u=p + h_u v_u$ with $h_u \in (h_1,h_2)$.
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
// p[2]
//    fixed point p=(g,G)
// v[2]
//    unstable vector
// lambda
//    unstable eigenvalue
// h
//    small increment in the direction of v
//    Keep it in sync with invmfld2!!
// a
//    axis $g=a$ parallel to the $G$ axis.
// piter
//    On exit, it contains the number of iterations of the poincare map
//    $\sixmap$ needed to take the unstable segment u_i to U_i and straddle
//    the $G$ axis.
// h_1,h_2
//    On exit, it contains the endpoints of the unstable segment u_i
//    bracketing the approximate root $p_u=p + h_u v_u$ with $h_u \in
//    (h_1,h_2)$.
// z[2]
//    On entry, it contains the previous intersection point (for the previous
//    energy level).
//    On exit, it contains the approximate intersection point for this energy
//    level.
// 
// RETURN VALUE
// ============
// Returns a non-zero error code to indicate an error and 0 to indicate
// success:
//    1: Problems computing the Poincare iterates.
//    2: No intersection found with the $G$ axis.
//
// NOTES
// =====
// For the moment, we work with the 3:1 resonant family of periodic orbits.
//
// If the first intersection point that we find is not part of the (continuous)
// family of primary intersections, then we keep looking for the primary one.
//
// Mental note: An alternative way to check that intersection belongs to
// primary family may be to check that piter is always increasing. If it was
// not, this would mean that the manifolds have developed a loop and a
// secondary intersection has appeared.
//
// We assume that fixed point p and st/unst manifolds are on the SEC2: {l=pi}
// section.
//
// Obs! parameter lambda_u is not used?

int 
approxint_del_unst (double mu, double H, int k, double p[2], double v[2],
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

   // 2. Discretize the linear segment between $p0=p+hv$ and $p1=P(p0)$ into
   // a set of NPOINTS points.

   // Compute $p_0$
   p0[0] = p[0] + h*v[0];
   p0[1] = p[1] + h*v[1];

   // Compute $p_1$
   p1[0] = p0[0];
   p1[1] = p0[1];
   status=prtbp_del_2d(mu,SEC2,H,k,p1,&ti);        // $p_1 = P(p_0)$
   if(status)
   {
      fprintf(stderr, "approxint_del: error computing Poincare map\n");
      return(1);
   }

   // Discretize linear segment
   disc(p0, p1, NPOINTS, l);

   // Make a copy of l in order to not modify it.
   memcpy(l_bak, l, 2*NPOINTS*sizeof(double));

   // Iterate the (discretized) linear segment "iter" times by the Poincare map
   for(iter=1;iter<=MAXITER;iter++)
   {
      status=u_i_del(mu, H, k, z, a, l, &i);
      if(status)
      {
	 fprintf(stderr, 
	       "approxint_del: error during %d-th iteration of linear segment\n",
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

   // update approximate intersection point
   z[0] = l[2*i];
   z[1] = l[2*i+1];
   return(0);
}

int 
approxint_del_st (double mu, double H, int k, double p[2], double v[2],
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

   // 2. Discretize the linear segment between $p0=p+hv$ and $p1=P(p0)$ into
   // a set of NPOINTS points.

   // Compute $p_0$
   p0[0] = p[0] + h*v[0];
   p0[1] = p[1] + h*v[1];

   // Compute $p_1$
   p1[0] = p0[0];
   p1[1] = p0[1];
   status=prtbp_del_2d_inv(mu,SEC2,H,k,p1,&ti);        // $p_1 = \sixmap^{-1}(p_0)$
   if(status)
   {
      fprintf(stderr, "approxint_del: error computing Poincare map\n");
      return(1);
   }

   // Discretize linear segment
   disc(p0, p1, NPOINTS, l);

   // Make a copy of l in order to not modify it.
   memcpy(l_bak, l, 2*NPOINTS*sizeof(double));

   // Iterate the (discretized) linear segment "iter" times by the Poincare map
   for(iter=1;iter<=MAXITER;iter++)
   {
      status=s_i_del(mu, H, k, z, a, l, &i);
      if(status)
      {
	 fprintf(stderr, 
	       "approxint_del: error during %d-th iteration of linear segment\n",
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

// name OF FUNCTION: u_i_del
//
// PURPOSE
// =======
// This function checks whether the n-th iteration does intersect the $G$
// axis or not. 
// If it does intersect the $G$ axis, it returns the unstable segment u_i
// such that its image $U_i$ crosses the $G$ axis.
//
// We discretize the unst domain into a set of NPOINTS segments u_1, u_2,
// ..., u_NPOINTS. 
// The image under iteration of this discretized version of the unst manifold
// is a set of segments U_1, U_2, ..., U_NPOINTS that approximates the
// nonlinear unst manifold.
// 
// We look for the first unst segment U_i that intersects the $G$ axis.
// Therefore, U_i contains an intersection point, and the segment u_i in the
// fundamental domain contains the preimage of an intersection point.
// It returns the unstable segment u_i such that its image $U_i$ crosses the
// $G$ axis.
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
//    axis $g=a$ parallel to the $G$ axis.
// lu
//    NPOINTS points approximating linear unstable fundamental domain.
//    On exit, this vector is modified with a new iteration of the
//    fundamental domain.
// idx
//    On exit, it contains the index $i$ corresponding to interval $u_i$
//    such that its image $U_i$ crosses the $G$ axis.
//    If the manifold does not cross the axis, we return i=-1.
// 
// RETURN VALUE
// ============
// Returns a non-zero error code to indicate an error and 0 to indicate
// success:
//    1: Problems computing the Poincare iterates.
//
// NOTES
// =====
// If the first intersection point that we find is not part of the (continuous)
// family of primary intersections, then we keep looking for the primary one.
//
// We assume that the st/unst manifolds are on the SEC2: {l=pi} section.

int 
u_i_del (double mu, double H, int k, double z[2], double a, double *l, int *idx)
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
	 if(prtbp_del_2d(mu,SEC2,H,k,l+2*i,&ti))
	 {
	    fprintf(stderr, "u_i_del: error computing Poincare map\n");
	    return(1);
	 }
   }

   // We look for the first unst segment U_i that crosses the $G$ axis.
   for(i=0; i<(NPOINTS-1); i++)
   {
      // Endpoints of segment U_i: 
      // P[2] = (g, G) = (l[2i], l[2i+1]),
      // Q[2] = (g, G) = (l[2i+2], l[2i+3]).
      if(((l[2*i]-a)*(l[2*i+2]-a)<=0) && 
	    fabs(l[2*i+1]-z[1]) < 0.02)	// point belongs to primary family
	 break;
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

int 
s_i_del (double mu, double H, int k, double z[2], double a, double *l, int *idx)
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
	 if(prtbp_del_2d_inv(mu,SEC2,H,k,l+2*i,&ti))
	 {
	    fprintf(stderr, "s_i_del: error computing Poincare map\n");
	    return(1);
	 }
   }

   // We look for the first st segment S_i that crosses the $G$ axis.
   for(i=0; i<(NPOINTS-1); i++)
   {
      // Endpoints of segment S_i: 
      // P[2] = (g, G) = (l[2i], l[2i+1]),
      // Q[2] = (g, G) = (l[2i+2], l[2i+3]).
      if(((l[2*i]-a)*(l[2*i+2]-a)<=0) && 
	    fabs(l[2*i+1]-z[1]) < 0.02)	// point belongs to primary family
	 break;
   }

   if(i==(NPOINTS-1))
   {
      // Manifold does not intersect $G$ axis
      *idx=-1;
      return(0);
   }

   // Manifold intersects the $G$ axis.

   //dx = l[2*i+3]-l[2*i+1];	// Q2 - P2
   //dy = l[2*i+2]-l[2*i];	// Q1 - P1
   //alpha2 = atan2(dy,dx);	
   //fprintf(stderr, "Approximate splitting half-angle: %f\n", alpha2);

   fprintf(stderr, "Approximate interval: %d\n", i);
   *idx=i;
   return(0);
}
