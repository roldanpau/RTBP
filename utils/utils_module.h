/*! 
  \file 
  \brief Utility functions.
  */

#include <stddef.h>	// size_t

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
