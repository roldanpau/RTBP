// ========================================
// Take points on section S1 to section S2.
// ========================================
// FILE:          $RCSfile: sec1sec2.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2013-02-20 10:02:48 $
//
// PURPOSE
// =======
// Take points on section S1 to section S2.
//
// NOTES
// =====
//
// OVERALL METHOD:
//
// 1. Input parameters from stdin:
// 
//    - mass parameter 
//
// 2. For each input line:
//
//    2.1. Input data:
//    - initial point "x" on section S1
//
//    2.2. Take x from section S1 to section S2 by the flow.
//
//    2.3. Output final point and integration time to stdout.

#include <stdio.h>
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE
#include <gsl/gsl_errno.h>	// gsl_set_error_handler_off
#include <rtbp.h>	// DIM
#include "sec1sec2_module.h"	// sec1sec2

int main( )
{
   double mu, ti;
   double x[DIM];
   int status;

   // Input mass parameter from stdin.
   if(scanf("%le", &mu) < 1)
   {
      perror("main: error reading input");
      exit(EXIT_FAILURE);
   }

   // Stop GSL default error handler from aborting the program
   gsl_set_error_handler_off();

   // For each input point x, do
   while(scanf("%le %le %le %le", x, x+1, x+2, x+3)==4)
   {
      // Take x from section S1 to section S2 by the flow.
      status=sec1sec2_inv(mu,x,&ti);
      if(status)
      {
	 
          fprintf(stderr, \
      "main: error taking point (%.15le,%.15le,%.15le,%.15le) from S1 to S2\n",
                  x[0],x[1],x[2],x[3]);
          exit(EXIT_FAILURE);
      }

      // Output final point and integration time to stdout.
      if(printf("%.15le %.15le %.15le %.15le %.15le\n", 
                  x[0], x[1], x[2], x[3], ti)<0)
      {
	 
          perror("main: error writting output");
          exit(EXIT_FAILURE);
      }
   }
   exit(EXIT_SUCCESS);
}
