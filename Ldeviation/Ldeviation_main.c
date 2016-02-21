/*! \file
    \brief Deviation of semi-major axis L wrt resonant one
    $Author: roldan $
    $Date: 2013-03-11 11:12:12 $
*/

#include <stdio.h>
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE
#include <gsl/gsl_errno.h>	// gsl_set_error_handler_off
#include <rtbp.h>	// DIM
#include "Ldeviation.h"

int main( )
{
   double mu;
   double H;		// energy level
   double x[DIM];
   double Ldev;		// maximum deviation
   int status;

   // Input mass parameter
   if(scanf("%le", &mu) < 1)
   {
      perror("main: error reading input");
      exit(EXIT_FAILURE);
   }

   // Stop GSL default error handler from aborting the program
   //gsl_set_error_handler_off();

   // Input initial conditions from stdin.
   while(scanf("%le %le %le %le %le", &H, x, x+1, x+2, x+3) == 5)
   {
      // Compute maximum deviation Ldev
      status=Ldeviation(mu,x,&Ldev);
      if(status)
      {
	 fprintf(stderr, \
	       "main: error computing maximum deviation\n");
	 exit(EXIT_FAILURE);
      }

      // Output H, Ldev to stdout.
      status = printf("%.15e %.15le\n", H, Ldev);
      if(status<0)
      {
	 perror("main: error writting output");
	 exit(EXIT_FAILURE);
      }
      fflush(NULL);
   }
   exit(EXIT_SUCCESS);
}
