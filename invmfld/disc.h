/*! \file
  \brief Discretize a linear segment into a given number of points
  */

/** 
  Discretize the linear segment between $p_0$ and $p_1$ into NPOINTS.

  \param[in] p0,p1
  endpoints of linear segment

  \param[in] n
  number of points in the discretization (n must be >=2)

  \param[out] l
  On return of the this function, it holds the discretization of the
  linear segment, i.e. an array of n 2D points.

  \returns
  a non-zero error code to indicate an error and 0 to indicate
  success.

  \pre
  Caller must make sure that enough space has been allocated for l (2*n
  doubles)

  \remark
  The points are equally spaced. 

  \remark
  The first point coincides with p0.
  The last point coincides with p1.
*/

int disc(double p0[2], double p1[2], int n, double *l);
