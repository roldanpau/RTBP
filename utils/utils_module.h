/*! 
  \file 
  \brief Utility functions.
  */

#include <stddef.h>	// size_t

extern const double TWOPI;

/** 
  Copy an array of doubles.

  \param[out] dst 	destination array
  \param[in] src 	source array
  \param[in] len	length of source array

  \return		destination array
  */
double * dblcpy(double * dst, double const * src, size_t len);

/** 
  Print an array of doubles.

  \param[in] x 		array
  \param[in] len	length of array

  \return		void
  */
void dblprint(double const *x, size_t len);

/** 
  Normalize angle in [-pi,pi).

  \param[in] fAng	angle

  \return		wrapped angle
  */
double WrapPosNegPI(double fAng);

/** 
  Normalize angle in [0,2pi).

  \param[in] fAng	angle

  \return		wrapped angle
  */
double WrapTwoPI(double fAng);
