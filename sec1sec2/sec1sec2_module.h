/*! \file
  \brief Take a point $p$ on section S1 to section S2 

   Consider the RTBP. 
   Let S1 be the Poincare section {y=0, vy>0}.
   Let S2 be the Poincare section {y=0, vy<0}.
   These functions take a point $p$ on section S1 to section S2 by either forward
   or backward integration.

  $Author: roldan $
  $Date: 2012-12-20 11:02:51 $
  */

#include <rtbp.h>   // DIM

/**
  This function takes a point $p$ on section S1 to section S2 by forward
  integration.

  \param[in] mu
      mass parameter for the RTBP
  \param[in,out] p
      Initial point, 4 coordinates: p=(x,y,p_x,p_y). 
      On return of the this function, it holds the image point on the section
      S2.
  \param[out] ti
      On exit, *ti holds the integration time to reach the section S2.

  \return 
  a non-zero error code to indicate an error:
  error computing poincare map,
  and 0 to indicate success.
  */

int sec1sec2(double mu, double p[DIM], double *ti);

/**
  This function takes a point $p$ on section S1 to section S2 by backward
  integration.

  \param[in] mu
      mass parameter for the RTBP
  \param[in,out] p
      Initial point, 4 coordinates: p=(x,y,p_x,p_y). 
      On return of the this function, it holds the image point on the section
      S2.
  \param[out] ti
      On exit, *ti holds the integration time to reach the section S2.

  \return 
  a non-zero error code to indicate an error:
  error computing poincare map,
  and 0 to indicate success.
  */

int sec1sec2_inv(double mu, double p[4], double *ti);
