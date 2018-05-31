/*! \file
    \brief Approximate Intersection of Invariant Manifolds

    $Author: roldan $
    $Date: 2013-03-26 22:10:03 $
*/

/** 
  Approximate intersection of unstable invariant manifold with symmetry line.

  Let $p$ be a hyperbolic fixed point for the 2D map 
  \f$\mathcal{P}: S \to S\f$, where S is the Poincare section SEC1 or SEC2, 
  corresponding to \f$\{l=0\}\f$ or \f$\{l=\pi\}\f$.

  Assume that \f$\lambda\f$ is the unstable eigenvalue, with \f$\lambda>1\f$.
  Let $v$ be the unstable eigenvector for the eigenvalue \f$\lambda\f$. 
  Let $W^u(p)$ be the unstable manifold of $p$.
  Let $g=a$ be a line parallel to the $G$ axis.

  This function computes an approximation to the first intersection of the
  unstable manifold with the line $g=a$ as we grow the manifold from the
  fixed point.

  We consider the unstable fundamental segment between the two points 
  $p+h_u v_u$ and $P(p+h_u v_u)$.
  We iterate the unstable fundamental domain until it intersects the line 
  $g=a$.
  
  We discretize the unst domain into a set of NPOINTS segments u_1, u_2,
  ..., u_NPOINTS. 
  The image under iteration of this discretized version of the unst manifold
  is a set of segments U_1, U_2, ..., U_NPOINTS that approximates the
  nonlinear unst manifold.
  We look for the first unst segment U_i that intersects the line $g=a$.
  Therefore, U_i contains an intersection point, and the segment u_i in the
  fundamental domains contain the preimage of an intersection point.

  This function returns the unstable segment u_i with endpoints $(h_1, h_2)$
  containing the approximate root $p_u=p + h_u v_u$ with \f$h_u \in
  (h_1,h_2)\f$.

  \remark
  We DO NOT assume that $p$ is located above the $p_x=0$ axis.

  \remark
  We DO NOT assume that $v=(x,p_x)$ points "to the right", i.e. 
  we DO NOT assume that the first component of $v$ is $x>0$.

  \param[in] mu         mass parameter for the RTBP
  \param[in] sec        Poincare section: sec={SEC1,SEC2}
  \param[in] H          energy value
  \param[in] k          number of iterates of the Poincare map
  \param[in] p          fixed point $p=(x,p_x)$
  \param[in] v          eigenvector associated to unstable direction
  \param[in] lambda     eigenvalue associated to unstable direction

  \param[in] h 		
  Small linear displacement along the manifold. This is typically obtained in
  function \ref h_opt.

  \param[in] br     branch type: br={LEFT, RIGHT}
  \param[in] a		line $g=a$ parallel to the $G$ axis.

  \param[out] piter
  On exit, it contains the number of iterations of the poincare map
  \f$\mathcal{P}\f$ needed to take the unstable segment u_i to U_i and
  straddle the line $g=a$.

  \param[out] h_1,h_2
  On exit, it contains the endpoints of the unstable segment u_i bracketing
  the approximate root $p_u=p + h_u v_u$ with \f$h_u \in (h_1,h_2)\f$.

  \param[out] z
  On exit, it contains the approximate intersection point for this energy
  level.

  \returns a non-zero error code to indicate an error and 0 to indicate
  success.

  \retval 1 Problems computing the Poincare iterates.
  \retval 2 No intersection found with the line $g=a$.
*/

int 
approxint_del_car_unst (double mu, section_t sec, double H, int k, 
        double p[2], double v[2], double lambda, double h, branch_t br, 
        double a, int *piter, double *h_1, double *h_2, double z[2]);

/**
  Approximate intersection of stable invariant manifold with symmetry line.

  Exactly as \ref approxint_del_car_unst.
  */
  
int
approxint_del_car_st (double mu, section_t sec, double H, int k, 
        double p[2], double v[2], double lambda, double h, branch_t br, 
        double a, int *piter, double *h_1, double *h_2, double z[2]);
