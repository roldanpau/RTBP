/*! \file
    \brief Poincare map of RTBP in Delaunay coordinates, integrating in Cartesian: main prog

    \author Pau Roldan
*/

#include <stdio.h>
#include <stdlib.h>	    // EXIT_SUCCESS, EXIT_FAILURE
#include <string.h>     // strcmp

#include <gsl/gsl_errno.h>	// gsl_set_error_handler_off

#include <rtbp.h>	        // DIM
#include "prtbp_del_car.h"	// section_t, prtbp_del_car

/**
  Poincare map of RTBP in Delaunay coordinates, but integrating the flow in Cartesian.

  Consider the RTBP in rotating Delaunay coordinates.
  Let S be the Poincare section SEC1, SEC2 or SECg, corresponding to {l=0},
  \f$\{l=\pi\}\f$, or \f$ \{g=0\} \f$.
  Let $x$ be a point in the section, and suppose that flow at $x$ is
  transversal to the section.
  Compute the n-th iterate of the Poincare map \f$P^n(x)\f$ of the RTBP.
  This procedure also computes, as a side product, the integration time to
  intersect the Poincare section "n" times.

  OVERALL METHOD:

  1. Input parameters from stdin:
  
    - mass parameter
    - type of section "sec"
    - number of iterates "n"

    A sequence of lines containing:

    - initial point in Delaunay "x_del"
    - initial point in Cartesian "x_car"
    
  2. Compute n-th iterate of the Poincare map, \f$P^n(x)\f$.

  3. Output final point $P^n(x)$ BOTH in Delaunay and Cartesian, 
  and integration time to stdout.
*/

int main( )
{
   double mu;
   section_t sec;
   double ti;		    // integration time
   double x_del[DIM];
   double x_car[DIM];
   int status, n;

   // auxiliary variables
   char section_str[10];        // holds input string "SEC1", "SEC2" etc

   // Input mass parameter, section, number of iterates from stdin.
   if(scanf("%le %s %d", &mu, section_str, &n) < 3)
   {
      perror("main: error reading input");
      exit(EXIT_FAILURE);
   }
   if (strcmp(section_str,"SEC1") == 0)
      sec = SEC1;
   else if (strcmp(section_str,"SEC2") == 0)
      sec = SEC2;
   else if (strcmp(section_str,"SECg") == 0)
      sec = SECg;
   else
   {
      perror("main: error reading section string");
      exit(EXIT_FAILURE);
   }

   // Stop GSL default error handler from aborting the program
   //gsl_set_error_handler_off();

   // Input initial conditions from stdin
   while( scanf("%le %le %le %le", x_del, x_del+1, x_del+2, x_del+3)==4 &&
           scanf("%le %le %le %le", x_car, x_car+1, x_car+2, x_car+3)==4 )
   {
       // Compute n-th iterate of Poincare map, $P^n(x)$.
       status=prtbp_del_car(mu,sec,n,x_del,x_car,&ti);
       if(status)
       {
          fprintf(stderr, \
            "main: error computing %d-th iterate of Poincare map\n",n);
          exit(EXIT_FAILURE);
       }

       // Output final point and integration time to stdout.
       status = printf("%.15le %.15le %.15le %.15le "
               "%.15le %.15le %.15le %.15le %.15le\n", \
         x_del[0], x_del[1], x_del[2], x_del[3], 
         x_car[0], x_car[1], x_car[2], x_car[3], ti);
       if(status<0)
       {
          perror("main: error writting output");
          exit(EXIT_FAILURE);
       }
   }
   exit(EXIT_SUCCESS);
}
