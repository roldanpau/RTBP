/*! 
  \file
  \brief Utility functions.
  */

#include <stdlib.h>	// malloc
#include <string.h>	// memcpy
#include <stdio.h>	// printf

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
