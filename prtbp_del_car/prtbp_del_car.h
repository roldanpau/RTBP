/*! \file
    \brief Poincare map of RTBP in Delaunay coordinates, but integrating the
    flow in Cartesian.
*/

#ifndef PRTBP_DEL_CAR_H_INCLUDED
#define PRTBP_DEL_CAR_H_INCLUDED

#include <section.h>	// section_t
#include "rtbp.h"	    // DIM

/// error bound (tolerance) for Poincare map
extern const double POINCARE_DEL_CAR_TOL;	

/**
  Poincare map of RTBP in Delaunay coordinates, but integrating the flow in
  Cartesian.

  Consider the RTBP in rotating Delaunay coordinates.
  Let S be the Poincare section SEC1, SEC2, or SECg, corresponding to {l=0},
  \f$\{l=\pi\}\f$, or \f$ \{g=0\} \f$.
  Let $x$ be a point in the section, and suppose that flow at $x$ is
  transversal to the section.
  Compute the n-th iterate of the Poincare map \f$P^n(x)\f$ of the RTBP.
  This procedure also computes, as a side product, the integration time to
  intersect the Poincare section "n" times.

  \param[in] mu mass parameter for the RTBP
  \param[in] sec type of Poincare section (sec = SEC1, SEC2, or SECg).

  \param[in] cuts 
  number of iterates of the Poincare map (cuts=n: n cuts with the Poincare
  section). This must be a positive (or zero) integer.

  \param[in,out] x 
  Initial point, 4 coordinates: (l, L, g, G). 
  On return of the this function, it holds the image point.

  \param[in,out] x_car
  Initial point in Cartesian, 4 coordinates: (x, y, p_x, p_y). 
  On return of the this function, it holds the image point.

  \param[out] ti
  On return, it holds the integration time to intersect the Poincare
  section "n" times.

  \return
  Returns a non-zero error code to indicate an error and 0 to indicate
  success.
  If an integration error is encountered, the function returns a non-zero
  value.

  \remark 
  The Poincare section is defined in Delaunay coordinates.
  However the flow is integrated in Cartesian. 

  \remark
  To avoid the change of coordinates from Delaunay to Cartesian, 
  which is numerically intensive (since one must solve Kepler's equation 
  \f$ u - e\sin u = l \f$), we pass the initial point BOTH in Delaunay 
  and Cartesian.
*/    

int prtbp_del_car(double mu, section_t sec, int cuts, double x[DIM], 
        double x_car[DIM], double *ti);

/**
  Inverse Poincare map of RTBP in Delaunay coordinates, 
  but integrating the flow in Cartesian.

  Consider the RTBP in rotating Delaunay coordinates.
  Let S be the Poincare section SEC1 or SEC2, corresponding to {l=0}
  or \f$\{l=\pi\}\f$.
  Let $x$ be a point in the section, and suppose that flow at $x$ is
  transversal to the section.
  Compute the n-th inverse iterate of the Poincare map \f$P^{-n}(x)\f$ 
  of the RTBP.
  This procedure also computes, as a side product, the integration time to
  intersect the Poincare section "n" times.

  \param[in] mu mass parameter for the RTBP
  \param[in] sec type of Poincare section (sec = SEC1 or SEC2).

  \param[in] cuts 
  number of iterates of the Poincare map (cuts=n: n cuts with the Poincare
  section). This must be a positive (or zero) integer.

  \param[in,out] x 
  Initial point, 4 coordinates: (l, L, g, G). 
  On return of the this function, it holds the image point.

  \param[in,out] x_car
  Initial point in Cartesian, 4 coordinates: (x, y, p_x, p_y). 
  On return of the this function, it holds the image point.

  \param[out] ti
  On return, it holds the integration time to intersect the Poincare
  section "n" times.
  Since we are computing the inverse Poincare map, this
  integration time "ti" must be negative.

  \return
  Returns a non-zero error code to indicate an error and 0 to indicate
  success.
  If an integration error is encountered, the function returns a non-zero
  value.

  \remark 
  The Poincare section is defined in Delaunay coordinates.
  However the flow is integrated in Cartesian. 

  \remark
  To avoid the change of coordinates from Delaunay to Cartesian, 
  which is numerically intensive (since one must solve Kepler's equation 
  \f$ u - e\sin u = l \f$), we pass the initial point BOTH in Delaunay 
  and Cartesian.
*/    

int prtbp_del_car_inv(double mu, section_t sec, int cuts, double x[DIM], 
        double x_car[DIM], double *ti);

#endif // PRTBP_DEL_CAR_H_INCLUDED
