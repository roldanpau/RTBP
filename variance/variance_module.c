/*! \file
  \brief Compute the first order of the variance $\sigma_0^2$.

  See the formula for the first order of the variance in Ansatz 2 of the second
  paper Stochastic3BP/Paper2NHIL3BP.

  */

#include <stdio.h>	// fprintf
#include <stdlib.h>	// EXIT_FAILURE
#include <complex.h>	// complex math type

double complex Mean(double complex C1, double complex C2)
{
	return (C1+C2)/2;
}
