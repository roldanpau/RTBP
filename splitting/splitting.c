/*! \file splitting.c
    \brief Splitting Angle of Invariant Manifolds

    $Author: roldan $
    $Date: 2013-03-11 11:36:04 $
*/

#include <stdio.h>
#include <math.h>	// sqrt
#include <rtbp.h>	// DIM

#include <utils_module.h>	// WrapPosNegPI
#include <hinv.h>
#include <prtbp_2d.h>
#include <prtbp_nl_2d_module.h>
#include <prtbp_nl.h>
#include <dprtbp_2d.h>	// dprtbp_nl_2d, dprtbp_2d_inv

int tanvec_u(double mu, double H, double v_u[2], int n, double p_u[2], 
      double w[2]);
int tanvec_s(double mu, double H, double v_s[2], int n, double p_s[2], 
      double w[2]);

/**
  Splitting Angle of Invariant Manifolds, computed using the unstable
  invariant manifold.

  Consider the 2D map \f$\mathcal{P}: \Sigma_- \to \Sigma_-\f$, which is
  assumed to be reversible with respect to the symmetry line $p_x=0$.
  Let $p$ be a hyperbolic fixed point for \f$\mathcal{P}\f$.
  For definiteness, we assume that $p$ is located above the $p_x=0$ axis.

  Assume that \f$\lambda\f$ is the unstable eigenvalue, with \f$\lambda>1\f$.
  Let $v$ be the unstable eigenvector for the eigenvalue \f$\lambda\f$. 
  For definiteness, we assume that $v=(x,p_x)$ points "to the right", i.e. we
  assume that the first component of $v$ is $x>0$, but we could use the other
  branch of the manifold.
  Let $W^u(p), W^s(p)$ be the unstable resp. stable manifold of $p$.

  Let $z=(x,0)$ be the "first" intersection point of the unstable manifold
  with the symmetry line $p_x=0$ as we grow the manifold from the fixed
  point.

  This function computes the splitting angle between the manifolds at the
  intersection point $z$.

  \param[in] mu         mass parameter for the RTBP
  \param[in] H          energy value
  \param[in] v          eigenvector associated to unstable direction

  \param[in] n
  number of desired iterations by the Poincare map in the unstable direction.

  \param[in] p[2]
  Preimage of the homoclinic point, i.e. point in the unstable fundamental
  segment such that
  \f$ z = \mathcal{P}^{n}(p) \f$.

  \param[out] angle
  On exit, it contains the splitting angle.

  \returns
  a non-zero error code to indicate an error and 0 to indicate success.
 */
//
// NOTES
// =====
// The splitting angle is computed in the following way: 
// Find the tangent vector, w_u, of the unstable manifold at z.
// Compute the half-angle $\alpha$ between w_u and the vertical vector (0,-1)
//      
// Finally, the angle between w_s and w_u is $2 \alpha$.

int splitting_unst(double mu, double H, double v[2], int n, double p[2],
      double *angle)
{
   double w[2]; // tangent vector to the unstable manifold at z
   //double dot;	  // dot product w_u\times w_s
   double alpha;  // splitting half-angle

   tanvec_u(mu, H, v, n, p, w);
   //printf("w_u: %.15le %.15le\n", w_u[0], w_u[1]);

   // Output splitting angle
   //dot = - w_u[1];

   alpha=atan2(-w[0],-w[1]);
   *angle = WrapPosNegPI(2.0*alpha);
   return(0);
}

/**
  Splitting Angle of Invariant Manifolds, computed using the stable
  invariant manifold.

  Exactly as \ref splitting_unst.
 */

int splitting_st(double mu, double H, double v[2], int n, double p[2],
      double *angle)
{
   double w[2]; // tangent vector to the stable manifold at z
   //double dot;	  // dot product w_u\times w_s
   double alpha;  // splitting half-angle

   tanvec_s(mu, H, v, n, p, w);
   //printf("w_u: %.15le %.15le\n", w_u[0], w_u[1]);

   // Output splitting angle
   //dot = - w_u[1];

   // For outer separatrix:
   alpha=atan2(-w[0],-w[1]);

   // For inner separatrix:
   //alpha=atan2(w_u[0],w_u[1]);
   //*angle = 2.0*acos(dot);
   *angle = 2.0*alpha;
   return(0);
}

// name OF FUNCTION: tanvec_u
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

int tanvec_u(double mu, double H, double v_u[2], int n, double p_u[2], 
      double w[2])
{
   double v[2]; 	// tangent vector to the manifold
   double x[DIM];	// point in the manifold
   double dp[4];	// Jacobian of 2d Poincare map
   double norm;		// norm of the tangent vector

   // auxiliary variables
   int i;
   double ti;

   v[0]=v_u[0]; 
   v[1]=v_u[1];
   x[0]=p_u[0]; x[1]=0; x[2]=p_u[1];
   // Compute x[3]=p_y by inverting the Hamiltonian.
   hinv(mu,SEC2,H,x);

   //printf("%le %le %le %le\n", x[0], x[1], v[0], v[1]);
   for(i=0; i<n; i++)
   {
      // x = P^i(p_u)
      // Jacobian of P at x
      dprtbp_nl_2d(mu,SEC2,4,x,dp);

      // w = DP*v
      w[0]=dp[0]*v[0]+dp[1]*v[1];
      w[1]=dp[2]*v[0]+dp[3]*v[1];

      // normalize w
      norm=sqrt(w[0]*w[0]+w[1]*w[1]);
      w[0]=w[0]/norm;
      w[1]=w[1]/norm;

      v[0]=w[0]; v[1]=w[1];

      prtbp_nl(mu,SEC2,4,x,&ti);
      //printf("%le %le %le %le\n", x[0], x[1], w[0], w[1]);
   }
   // On exit, we have:
   //    x = z = P^{i+1}(p_u)
   //    w is the tangent vector of the manifold at z
   return(0);
}

int tanvec_s(double mu, double H, double v_s[2], int n, double p_s[2], 
      double w[2])
{
   double v[2]; 	// tangent vector to the manifold
   double x[2];		// point in the manifold
   double dp[4];	// Jacobian of 2d Poincare map
   double norm;		// norm of the tangent vector

   // auxiliary variables
   int i;
   double ti;

   v[0]=v_s[0]; 
   v[1]=v_s[1];
   x[0]=p_s[0]; x[1]=p_s[1];
   //printf("%le %le %le %le\n", x[0], x[1], v[0], v[1]);
   for(i=0; i<n; i++)
   {
      // x = P^{-i}(p_s)
      // Jacobian of P at x
      dprtbp_2d_inv(mu,SEC2,H,2,x,dp);

      // w = DP*v
      w[0]=dp[0]*v[0]+dp[1]*v[1];
      w[1]=dp[2]*v[0]+dp[3]*v[1];

      // normalize w
      norm=sqrt(w[0]*w[0]+w[1]*w[1]);
      w[0]=w[0]/norm;
      w[1]=w[1]/norm;

      v[0]=w[0]; v[1]=w[1];

      prtbp_2d_inv(mu,SEC2,H,2,x,&ti);
      //printf("%le %le %le %le\n", x[0], x[1], w[0], w[1]);
   }
   // On exit, we have:
   //    x = z = P^{-(i+1)}(p_s)
   //    w is the tangent vector of the manifold at z
   return(0);
}

