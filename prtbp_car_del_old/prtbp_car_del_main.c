/*! \file
    \brief Poincare map of RTBP in Delaunay coordinates, integrating in Cartesian: main prog

    \author Pau Roldan
*/

#include <stdio.h>
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE
#include <string.h>     // strcmp

#include <gsl/gsl_errno.h>	// gsl_set_error_handler_off

#include <rtbp.h>	// DIM
#include "prtbp_car_del.h"	// section_t, prtbp_car_del

/**
  Poincare map of RTBP in Delaunay coordinates, integrating in Cartesian: main prog

  Consider the RTBP in rotating coordinates.
  Let S be the Poincare section SEC1 or SEC2, corresponding to {l=0}
  or \f$\{l=\pi\}\f$.
  Let $x$ be a point in the section, and suppose that flow at $x$ is
  transversal to the section.
  Compute the n-th iterate of the Poincare map $P^n(x)$ of the RTBP.
  This procedure also computes, as a side product, the integration time to
  intersect the Poincare section "n" times.

  OVERALL METHOD:

  1. Input parameters from stdin:
  
    - mass parameter
    - type of section "sec"
    - number of iterates "n"
    - initial point "x"
    
  2. Compute n-th iterate of the Poincare map, \f$P^n(x)\f$.

  3. Output final point $P^n(x)$ and integration time to stdout.
*/

int main( )
{
   double mu;
   section_t sec;
   double ti;		// integration time
   double x[DIM];
   int status, n;

   // auxiliary variables
   char section_str[10];        // holds input string "SEC1", "SEC2" etc

   // Input mass parameter, section, number of iterates and initial condition from stdin.
   if(scanf("%le %s %d %le %le %le %le", &mu, section_str, &n, x, x+1, x+2, x+3) < 7)
   {
      perror("main: error reading input");
      exit(EXIT_FAILURE);
   }

   if (strcmp(section_str,"SEC1") == 0)
      sec = SEC1;
   else if (strcmp(section_str,"SEC2") == 0)
      sec = SEC2;
   else
   {
      perror("main: error reading section string");
      exit(EXIT_FAILURE);
   }

   // Stop GSL default error handler from aborting the program
   //gsl_set_error_handler_off();

   // Compute n-th iterate of Poincare map, $P^n(x)$.
   status=prtbp_car_del(mu,sec,n,x,&ti);
   if(status)
   {
      fprintf(stderr, \
	    "main: error computing %d-th iterate of Poincare map\n",n);
      exit(EXIT_FAILURE);
   }

   // Output final point and integration time to stdout.
   status = printf("%.15le %.15le %.15le %.15le %.15le\n", \
	 x[0], x[1], x[2], x[3], ti);
   if(status<0)
   {
      perror("main: error writting output");
      exit(EXIT_FAILURE);
   }
   exit(EXIT_SUCCESS);
}
