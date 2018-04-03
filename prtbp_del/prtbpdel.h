/*! \file prtbpdel.h
    \brief Poincare map of RTBP in Delaunay coordinates

    $Author: roldan $
    $Date: 2013-03-26 22:23:08 $
*/

#ifndef PRTBPDEL_H_INCLUDED
#define PRTBPDEL_H_INCLUDED

#include <stdbool.h>	// bool
#include <rtbp.h>	    // DIM
#include <section.h>	// section_t

bool onsection_del (section_t sec, double x[DIM]);
bool crossing_fwd_del (section_t sec, double x[DIM], double y[DIM]);
bool crossing_bwd_del (section_t sec, double x[DIM], double y[DIM]);

/**
  Poincare map of RTBP in Delaunay coordinates

  Consider the RTBP in rotating coordinates.
  Let S be the Poincare section SEC1 or SEC2, corresponding to {l=0} or
  \f$\{l=\pi\}\f$.
  Let $x$ be a point in the section, and suppose that flow at a point $x$ is 
  transversal to the section.
  Compute the n-th iterate of the Poincare map $P^n(x)$ of the RTBP.
  This procedure also computes, as a side product, the integration time to
  intersect the Poincare section "n" times.

  \param[in] mu mass parameter for the RTBP
  \param[in] sec type of Poincare section (sec = SEC1 or SEC2).

  \param[in] cuts 
  number of iterates of the Poincare map (cuts=n: n cuts with the Poincare
  section). This must be a positive (or zero) integer.

  \param[in,out] x 
  Initial point, 4 coordinates: (l, L, g, G). 
  On return of the this function, it holds the image point $P^n(x)$.

  \param[out] ti
  On return, it holds the integration time to intersect the Poincare
  section "n" times.

  \return
  Returns a non-zero error code to indicate an error and 0 to indicate
  success.
  If an integration error is encountered, the function returns a non-zero
  value.

  \remark 
  To compute the inverse Poincare map (i.e., negative iterates of the
  Poincare map), see \ref prtbp_del_inv.

  \remark
  On successful return of this function, the point $x$ is exactly on the
  section, i.e. we set coordinate $l$ exactly equal to zero/pi.

  \remark
  Since $g$ is an angle, we normalize it to \f$[0,2\pi)\f$.
  */

int prtbp_del(double mu, section_t sec, int cuts, double x[DIM], double *ti);

/**
  Inverse Poincare map of RTBP in Delaunay coordinates

  Consider the RTBP in rotating coordinates.
  Let S be the Poincare section SEC1 or SEC2, corresponding to {l=0} or
  \f$\{l=\pi\}\f$.
  Let $x$ be a point in the section, and suppose that flow at a point $x$ is 
  transversal to the section.
  Compute the n-th iterate of the inverse Poincare map $P^{-n}(x)$ of the RTBP.
  This procedure also computes, as a side product, the integration time to
  intersect the Poincare section "n" times.

  \param[in] mu mass parameter for the RTBP
  \param[in] sec type of Poincare section (sec = SEC1 or SEC2).

  \param[in] cuts 
  number of iterates of the Poincare map (cuts=n: n cuts with the Poincare
  section). This must be a POSITIVE (or zero) integer.

  \param[in,out] x 
  Initial point, 4 coordinates: (l, L, g, G). 
  On return of the this function, it holds the image point $P^{-n}(x)$.

  \param[out] ti
  On return, it holds the integration time to intersect the Poincare
  section "n" times.

  \return
  Returns a non-zero error code to indicate an error and 0 to indicate
  success.
  If an integration error is encountered, the function returns a non-zero
  value.

  \remark 
  To compute the direct Poincare map (i.e., positive iterates of the
  Poincare map), see \ref prtbp_del.

  \remark
  On successful return of this function, the point $x$ is exactly on the
  section, i.e. we set coordinate $l$ exactly equal to zero/pi.

  \remark
  Since $g$ is an angle, we normalize it to \f$[0,2\pi)\f$.
  */

int prtbp_del_inv(double mu, section_t sec, int cuts, double x[DIM], double *ti);

#endif // PRTBPDEL_H_INCLUDED
