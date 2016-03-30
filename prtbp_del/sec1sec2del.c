// ============================================
// Take points on section SEC1 to section SEC2.
// ============================================
// FILE:          $RCSfile: sec1sec2del_main.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2012-11-26 10:42:05 $
//
// PURPOSE
// =======
// Take points on section SEC1 to section SEC2.
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
//    - eccentricity "e"
//    - energy value "H"
//    - initial point "x" on section SEC1
//
//    2.2. Take x from section SEC1 to section SEC2 by the flow.
//
//    2.3. Output final point and integration time to stdout.

#include <stdio.h>
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE
#include <gsl/gsl_errno.h>	// gsl_set_error_handler_off
#include "prtbpdel.h"		// SEC1, SEC2
#include "prtbpdel_2d.h"	// sec1sec2_del

int main( )
{
   double mu, e, H, ti;
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

   // For each eccentricity e in the range, do
   while(scanf("%le %le %le %le", &e, &H, x, x+1)==4)
   {
      // Take x from section SEC1 to section SEC2 by the flow.
      status=sec1sec2_del(mu,H,x,&ti);
      if(status)
      {
	 fprintf(stderr, \
	       "main: error taking point (%.15e,%.15e) from SEC1 to SEC2\n",
	       x[0],x[1]);
	 exit(EXIT_FAILURE);
      }

      // Output final point and integration time to stdout.
      if(printf("%e %.15le %.15le %.15le %.15le\n", e, H, x[0], x[1], ti)<0)
      {
	 perror("main: error writting output");
	 exit(EXIT_FAILURE);
      }
   }
   exit(EXIT_SUCCESS);
}
