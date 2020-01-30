/*! 
  \file
  \brief Utility functions.
  */

#include <stdlib.h>	// malloc
#include <string.h>	// memcpy
#include <stdio.h>	// printf
#include <math.h>	// M_PI, floor

const double TWOPI = 2*M_PI;

double *dblcpy(double * dst, double const * src, size_t len)
{
   memcpy(dst, src, len * sizeof(double));
   return dst;
}

void dblprint(double const *x, size_t len)
{
   int i;
   for(i=0; i<len; i++)
      printf("%.15le ", x[i]);
}

// Floating-point modulo
// The result (the remainder) has same sign as the divisor.
// Similar to matlab's mod(); 
// Not similar to fmod() -   Mod(-3,4)= 1   fmod(-3,4)= -3
double Mod(double x, double y)
{
    if(y == 0)
        return x;

    double m= x - y * floor(x/y);
    return m;
}

// wrap [rad] angle to [-pi,pi)
double WrapPosNegPI(double fAng)
{
    return Mod(fAng + M_PI, TWOPI) - M_PI;
}

// wrap [rad] angle to [0..2\pi)
double WrapTwoPI(double fAng)
{
    return Mod(fAng, TWOPI);
}

// L_2 norm of an array x of size n
double l2_norm(double const* x, int n)
{
	double accum = 0.;
	for (int i = 0; i<n; ++i)
	{
		accum += x[i] * x[i];
	}
	return sqrt(accum);
}

