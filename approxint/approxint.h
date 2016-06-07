/*! \file
    \brief Approximate Intersection of Invariant Manifolds

    $Author: roldan $
    $Date: 2013-03-26 22:10:03 $
*/

/** 
  Approximate intersection of unstable invariant manifold with symmetry line.

  Let $p$ be a hyperbolic fixed point for the 2D map \f$\mathcal{P}: \Sigma_-
  \to \Sigma_-\f$.
  For definiteness, we assume that $p$ is located above the $p_x=0$ axis.

  Assume that \f$\lambda\f$ is the unstable eigenvalue, with \f$\lambda>1\f$.
  Let $v$ be the unstable eigenvector for the eigenvalue \f$\lambda\f$. 
  For definiteness, we assume that $v=(x,p_x)$ points "to the right", i.e. we
  assume that the first component of $v$ is $x>0$, but we could use the other
  branch of the manifold.
  Let $W^u(p)$ be the unstable manifold of $p$.
  Let $p_x=a$ be a line parallel to the $x$ axis.

  This function computes an approximation to the "first" intersection of the
  unstable manifold with the line $p_x=a$ as we grow the manifold from the
  fixed point.

  We consider the unstable fundamental segment between the two points 
  $p+h_u v_u$ and $P(p+h_u v_u)$.
  We iterate the unstable fundamental domain until it intersects the $x$
  axis.
  
  We discretize the unst domain into a set of NPOINTS segments u_1, u_2,
  ..., u_NPOINTS. 
  The image under iteration of this discretized version of the unst manifold
  is a set of segments U_1, U_2, ..., U_NPOINTS that approximates the
  nonlinear unst manifold.
  We look for the first unst segment U_i that intersects the $x$ axis.
  Therefore, U_i contains an intersection point, and the segment u_i in the
  fundamental domains contain the preimage of an intersection point.

  This function returns the unstable segment u_i with endpoints $(h_1, h_2)$
  containing the approximate root $p_u=p + h_u v_u$ with \f$h_u \in
  (h_1,h_2)\f$.

  \remark If the first intersection point that we find is not part of the
  (continuous) family of primary intersections, then we keep looking for the
  primary one.

  \param[in] mu         mass parameter for the RTBP
  \param[in] H          energy value
  \param[in] k          number of iterates of the Poincare map
  \param[in] p          fixed point $p=(x,p_x)$
  \param[in] v          eigenvector associated to unstable direction
  \param[in] lambda     eigenvalue associated to unstable direction

  \param[in] h 		
  Small linear displacement along the manifold. This is typically obtained in
  function \ref h_opt.

  \param[in] a		line $p_x=a$ parallel to the $x$ axis.

  \param[out] piter
  On exit, it contains the number of iterations of the poincare map
  \f$\mathcal{P}\f$ needed to take the unstable segment u_i to U_i and
  straddle the $x$ axis.

  \param[out] h_1,h_2
  On exit, it contains the endpoints of the unstable segment u_i bracketing
  the approximate root $p_u=p + h_u v_u$ with \f$h_u \in (h_1,h_2)\f$.

  \param[in,out] z
  On entry, it contains the previous intersection point (for the previous
  energy level).
  On exit, it contains the approximate intersection point for this energy
  level.

  \returns a non-zero error code to indicate an error and 0 to indicate
  success.

  \retval 1 Problems computing the Poincare iterates.
  \retval 2 No intersection found with the $x$ axis.
*/
int 
approxint_unst (double mu, double H, int k, double p[2], double v[2],
      double lambda, double h, double a, 
      int *piter, double *h_1, double *h_2, double z[2]);

/**
  Approximate intersection of stable invariant manifold with symmetry line.

  Exactly as \ref approxint_unst.
  */
  
int 
approxint_st (double mu, double H, int k, double p[2], double v[2],
      double lambda, double h, double a, 
      int *piter, double *h_1, double *h_2, double z[2]);
