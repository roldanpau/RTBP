/*! \file cardel_2d.h
    \brief Cartesian to Delaunay Change of Coordinates (2D)

    Obtain rotating Delaunay coordinates (l,L,g,G) from rotating cartesian
    (x,px).

    $Author: roldan $
    $Date: 2013-03-26 22:15:14 $
*/

#include <rtbp.h>	//DIM

/**
  Obtain rotating Delaunay coordinates (l,L,g,G) from rotating cartesian
  (x,px).

  \remark
  WE ASSUME that input point is on the Poincare section 
  \f$ \Sigma_- = {y=0, \dot y<0} \f$. We can recover the 4 euclidean coordinates
  (x,y,px,py) using the energy relation.

  \param[in] mu	mass parameter for the RTBP
  \param[in] H 	energy value
  \param[in] z 	point in rotating cartesian coordinates, z=(x,p_x).
  \param[out] Y point in rotating Delaunay coordinates, Y=(l,L,g,G).

  \return a non-zero error code to indicate an error and 0 to indicate
  success.
*/    
int cardel_2d(double mu, double H, double z[2], double Y[DIM]);
