/*! \file
    \brief Hyperbolic Splitting (Eigenvalues/vects) Associated to a Fixed
    Point

    $Author: roldan $
    $Date: 2013-03-11 11:22:18 $
*/

#include <prtbp.h>	// section_t

/** Error computing derivative of 2D Poincare map of RTBP. */
extern const int ERR_DPRTBP_2D;

/**  Eigenvalues are not real, so fixed point is not hyperbolic! */
extern const int ERR_NOTREAL;

/**
  Hyperbolic Splitting (Eigenvalues/vects) Associated to a Fixed
  Point.

  Consider a 2D map $F$ and its derivative $DF$.
  Given a hyperbolic fixed point $p$ of $F$, compute its hyperbolic
  splitting, i.e. its unstable/stable eigenvalues and eigenvectors. 
  Since $DF$ is a 2x2 matrix, we compute this simply by hand, see for
  instance
  http://www.math.harvard.edu/archive/21b_fall_04/exhibits/2dmatrices/index.html

  \remark
  Eigenvectors are normalised to unit magnitude.

  \remark 
  Eigenvectors are chosen to point "to the right" of the fixed point, i.e.
  $v=(x,p_x)$ with first component $x>0$.

  \param[in] mu         mass parameter of RTBP.
  \param[in] sec        Poincare section
  \param[in] H          energy value of periodic orbit associated to $p$
  \param[in] n          number of cuts with section

  \param[in] p          
  Argument of the 2d derivative, 2 coordinates: p=(x,p_x).

  \param[out] eval
  On exit, eval[0] holds the unstable eigenvalue, and eval[1] holds the
  stable eigenvalue.

  \param[out] evec
  On exit, the first row evec[0,1] holds the unstable eigenvector, and the
  second row evect[2,3] holds the stable eigenvector.

  \returns a non-zero error code to indicate an error and 0 to indicate
  success.
 
  \retval ERR_DPRTBP_2D
  Error computing derivative of 2D Poincare map of RTBP.
  
  \retval ERR_NOTREAL
  Eigenvalues are not real, so fixed point is not hyperbolic!
*/

int hyper(double mu, section_t sec, double H, int n, double p[2], double eval[2], 
      double evec[4]);
