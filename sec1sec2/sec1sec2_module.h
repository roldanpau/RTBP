/*! \file
  \brief Take a point $p$ on section S1 to section S2 

   Consider the RTBP. 
   Let S1 be the Poincare section {y=0, vy>0}.
   Let S2 be the Poincare section {y=0, vy<0}.
   Fix the Hamiltonian to a given value "H". 
   Using this energy condition, we can work with only two variables, $(x,p_x)$.
   The third variable $y$ is 0 since we look at the Poincare section, and the
   fourth variable $p_y$ can be obtained from the energy condition.
   These functions take a point $p$ on section S1 to section S2 by either forward
   or backward integration.

  $Author: roldan $
  $Date: 2012-12-20 11:02:51 $
  */

/**
  This function takes a point $p$ on section S1 to section S2 by forward
  integration.

  \param[in] mu
      mass parameter for the RTBP
  \param[in] H
      energy value
  \param[in,out] p
      Initial point, 2 coordinates: p=(x,p_x). 
      On return of the this function, it holds the image point on the section
      S2.
  \param[out] ti
      On exit, *ti holds the integration time to reach the section S2.

  \return 
  a non-zero error code to indicate an error:
  error inverting the Hamiltonian, 
  error computing poincare map,
  and 0 to indicate success.
  */

int sec1sec2(double mu, double H, double p[2], double *ti);

/**
  This function takes a point $p$ on section S1 to section S2 by backward
  integration.

  \param[in] mu
      mass parameter for the RTBP
  \param[in] H
      energy value
  \param[in,out] p
      Initial point, 2 coordinates: p=(x,p_x). 
      On return of the this function, it holds the image point on the section
      S2.
  \param[out] ti
      On exit, *ti holds the integration time to reach the section S2.

  \return 
  a non-zero error code to indicate an error:
  error inverting the Hamiltonian, 
  error computing poincare map,
  and 0 to indicate success.
  */

int sec1sec2_inv(double mu, double H, double p[2], double *ti);
