// ========================================
// Take points on section S2 to section S1.
// ========================================
// FILE:          $RCSfile: sec2sec1.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2013-02-20 10:02:48 $
//
// PURPOSE
// =======
// Take points on section S2 to section S1.
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
//    - energy value "H"
//    - initial point "x" on section S2
//
//    2.2. Take x from section S2 to section S1 by the flow.
//
//    2.3. Output final point and integration time to stdout.

#include <stdio.h>
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE
#include <gsl/gsl_errno.h>	// gsl_set_error_handler_off
#include "sec2sec1_module.h"	// sec2sec1

int main( )
{
   double mu, H, ti;
   double x[2];
   int status;

   // Input mass parameter from stdin.
   if(scanf("%le", &mu) < 1)
   {
      perror("main: error reading input");
      exit(EXIT_FAILURE);
   }

   // Stop GSL default error handler from aborting the program
   gsl_set_error_handler_off();

   // For each energy level H in the range, do
   while(scanf("%le %le %le", &H, x, x+1)==3)
   {
      // Take x from section S2 to section S1 by the flow.
      status=sec2sec1(mu,H,x,&ti);
      if(status)
      {
	 fprintf(stderr, \
	       "main: error taking point (%.15e,%.15e) from S2 to S1\n",
	       x[0],x[1]);
	 exit(EXIT_FAILURE);
      }

      // Output final point and integration time to stdout.
      if(printf("%.15e %.15le %.15le %.15le\n", H, x[0], x[1], ti)<0)
      {
	 perror("main: error writting output");
	 exit(EXIT_FAILURE);
      }
   }
   exit(EXIT_SUCCESS);
}
