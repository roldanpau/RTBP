// ======================================
// Splitting Angle of Invariant Manifolds
// ======================================
// FILE:          $RCSfile: splittingdel.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2012-06-06 13:07:51 $
//
// FUNCTIONS
// =========
//
// splitting_del
// -------------
// Consider the Restricted Three Body Problem. Assume that energy is fixed:
// $H(l,L,g,G)=\bar H$.
// Let SEC2: {l=pi} be a Poincare section, and $P(g,G)$ be the associated 2D
// Poincare map.
// Let $p(g,G)$ be a fixed point of the Poincare map associated to a
// periodic orbit of the flow.
// Let $W^u(p), W^s(p)$ be the unstable/stable manifold of $p$.
// Let $z$ be the intersection point $z\in W^u(p) \cap W^s(p)$ that is closer
// to the "fold" of the manifolds. 
// This function computes the splitting angle between the manifolds at the
// intersection point $z$.
//
// tanvec_u_del
// --------
// This function computes the tangent vector to the unstable manifold at the
// homoclinic point z.


#include <stdio.h>
#include <math.h>	// sqrt

#include <prtbpdel_2d.h>	// prtbp_del_2d, prtbp_del_2d_inv
#include <dprtbpdel_2d.h>	// dprtbp_del_2d, dprtbp_del_2d_inv

int tanvec_u_del(double mu, double H, double v_u[2], int n, double p_u[2], 
      double w[2]);

// name OF FUNCTION: splitting_del
//
// PURPOSE
// =======
// Consider the Restricted Three Body Problem. Assume that energy is fixed:
// $H(l,L,g,G)=\bar H$.
// Let SEC2: {l=pi} be a Poincare section, and $P(g,G)$ be the associated 2D
// Poincare map.
// Let $p(g,G)$ be a fixed point of the Poincare map associated to a
// periodic orbit of the flow.
// Let $W^u(p), W^s(p)$ be the unstable/stable manifold of $p$.
// Let $z$ be the intersection point $z\in W^u(p) \cap W^s(p)$ that is closer
// to the "fold" of the manifolds. 
// This function computes the splitting angle between the manifolds at the
// intersection point $z$.
//
// The splitting angle is computed in the following way: 
// Find the tangent vector, w_u, of the unstable manifold at z.
// Compute the half-angle $\alpha$ between w_u and the horizontal vector (1,0)
//      
// Finally, the angle between w_s and w_u is $2 \alpha$.
//
// PARAMETERS
// ==========
// mu
//    mass parameter for the RTBP
// H
//    energy value
// v_u[2]
//    unstable vector
// n
//    number of desired iterations by the Poincare map in the unstable
//    direction.
// p_u[2]
//    Points in the unstable segments such that
//       P^{n}(p_u).
// w_u[2]
//    On exit, it contains tangent vector to the manifold at homocl. pt.
// angle
//    On exit, it contains the splitting angle.
// 
// RETURN VALUE
// ============
// Returns a non-zero error code to indicate an error and 0 to indicate
// success:
//
// NOTES
// =====

int splitting_del(double mu, double H, double v_u[2], int n, double p_u[2],
      double w_u[2], double *angle)
{
   double dot;	  // dot product w_u\times w_s
   double alpha;  // splitting half-angle

   tanvec_u_del(mu, H, v_u, n, p_u, w_u);

   // Output splitting angle
   //dot = - w_u[1];

   // For outer separatrix:
   //alpha=atan2(-w_u[0],-w_u[1]);

   // For inner separatrix:
   alpha=atan2(w_u[1],w_u[0]);
   //*angle = 2.0*acos(dot);
   *angle = 2.0*alpha;
   return(0);
}

// name OF FUNCTION: tanvec_u_del
//
// PURPOSE
// =======
// This function computes the tangent vector to the unstable manifold at the
// homoclinic point z.
//
// Let $p_u = p+h_u v_u$ be the point in the unstable segment, such that
//    $z=P^n(p_u)$.
//
// Consider the tangent vector v_u to the manifold at p_u.  (Recall that at
// this point the linear approximation is good enough, so we can take as v_u
// the unstable eigenvector at the fixed point p.) 
// Multiply v_u by the Jacobian of P at the successive iterates
// P^i(p_u), for i=0,...,n-1.
// This way, we obtain the tangent vector, say w, of the manifold at z.
//
// PARAMETERS
// ==========
// mu
//    mass parameter for the RTBP
// H
//    energy value
// v_u
//    unstable vector
// n
//    number of desired iterations by the Poincare map in the unstable
//    direction.
// p_u[2]	
//    Point in the unstable segment such that
//       z = P^n(p_u).
// w
//    On exit, it contains the tangent vector of the manifold at z.
// 
// RETURN VALUE
// ============
// Returns a non-zero error code to indicate an error and 0 to indicate
// success.
//
// NOTES
// =====
// Since the tangent vector is expanded a lot due to multiplication by the
// Jacobian, we normalize it to norm 1: ||w||=1
//
// We assume that we are interested on the 3:1 resonance.

int tanvec_u_del(double mu, double H, double v_u[2], int n, double p_u[2], 
      double w[2])
{
   double v[2]; 	// tangent vector to the manifold
   double x[2];		// point in the manifold
   double dp[4];	// Jacobian of 2d Poincare map
   double norm;		// norm of the tangent vector

   // auxiliary variables
   int i;
   double ti;

   v[0]=v_u[0]; 
   v[1]=v_u[1];
   x[0]=p_u[0]; x[1]=p_u[1];
   //printf("%le %le %le %le\n", x[0], x[1], v[0], v[1]);
   for(i=0; i<n; i++)
   {
      // x = P^i(p_u)
      // Jacobian of P at x
      dprtbp_del_2d(mu,SEC2,H,3,x,dp);

      // w = DP*v
      w[0]=dp[0]*v[0]+dp[1]*v[1];
      w[1]=dp[2]*v[0]+dp[3]*v[1];

      // normalize w
      norm=sqrt(w[0]*w[0]+w[1]*w[1]);
      w[0]=w[0]/norm;
      w[1]=w[1]/norm;

      v[0]=w[0]; v[1]=w[1];

      prtbp_del_2d(mu,SEC2,H,3,x,&ti);
      //printf("%le %le %le %le\n", x[0], x[1], w[0], w[1]);
   }
   // On exit, we have:
   //    x = z = P^{i+1}(p_u)
   //    w is the tangent vector of the manifold at z
   return(0);
}

