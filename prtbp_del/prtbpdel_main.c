/*! \file prtbpdel_main.c
    \brief Poincare map of RTBP in Delaunay coordinates: main prog.

    $Author: roldan $
    $Date: 2013-03-26 22:23:08 $
*/

#include <stdio.h>
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE
#include <string.h>	// strcmp
#include <gsl/gsl_errno.h>	// gsl_set_error_handler_off
#include <rtbp.h>	// DIM
#include <section.h>
#include "prtbpdel.h"

/**
  Poincare map of RTBP in Delaunay coordinates: main prog.

  Let "sec" be a Poincare section in the RTBP, where
  	- sec = SEC1 means section {l=0}
  	- sec = SEC2 means section {l=pi}.
  	- sec = SECg means section {g=0}.
  	- sec = SECg2 means section {g=pi}.
  This program computes the n-th iterate of the Poincare map, $P^n(x)$.
  NOTE: n<0 means use inverse Poincare map!
 
  1. Input params (stdin):

  - mass parameter
  - type of section "sec"
  - number of iterates "n". 

  2. for each initial point "x" from stdin, do
     2.1. Compute n-th iterate of the Poincare map, $P^n(x)$.
     2.2. Output final point $P^n(x)$ and integration time to stdout.

  \remark
  We do not check (and thus we do not impose) that the initial point "x" is
  on the Poincare section.
 
 */

int main( )
{
   double mu;
   section_t sec;
   double ti;		// integration time
   double x[DIM];
   int status, n;

   // auxiliary variables
   char section_str[10];	// holds input string "SEC1", "SEC2" etc

   // Input mass parameter, section type, number of iterates from stdin.
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
   else if (strcmp(section_str,"SECg2") == 0)
      sec = SECg2;
   else 
   {
      perror("main: error reading section string");
      exit(EXIT_FAILURE);
   }

   // Stop GSL default error handler from aborting the program
   //gsl_set_error_handler_off();

   // Input initial conditions from stdin.
   while(scanf("%le %le %le %le", x, x+1, x+2, x+3) == 4)
   {
	   if(n>0)
	   {
		  // Compute n-th iterate of Poincare map, $P^n(x)$.
		  status=prtbp_del(mu,sec,n,x,&ti);
		  if(status)
		  {
		 fprintf(stderr, \
			   "main: error computing %d-th iterate of Poincare map\n",n);
		 exit(EXIT_FAILURE);
		  }
	   }
	   else if(n<0)
	   {
		  // Compute n-th iterate of inverse Poincare map, $P^{-n}(x)$.
		  status=prtbp_del_inv(mu,sec,-n,x,&ti);
		  if(status)
		  {
		 fprintf(stderr, \
			   "main: error computing %d-th inverse iterate of Poincare map\n",
			   -n); 
		 exit(EXIT_FAILURE);
		  }
	   }

      // Output final point and integration time to stdout.
      status = printf("%.15le %.15le %.15le %.15le %.15le\n", \
	    x[0], x[1], x[2], x[3], ti);
      if(status<0)
      {
	 perror("main: error writting output");
	 exit(EXIT_FAILURE);
      }
   }
   exit(EXIT_SUCCESS);
}
