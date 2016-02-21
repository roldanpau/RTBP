// ==========================================
// Continuation of the shift A^+ wrt energy H
// ==========================================
// FILE:          $RCSfile: inner_ell_main.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2011-07-01 09:20:28 $
//
// PURPOSE:
// Consider the family of 7:1 resonant periodic orbits in the RTBP.
// Compute the periodic orbit for a given (wide) range of energies, and
// compute how the integral $A^+$ changes when we change the energy H.
// This function determines the twist of the inner map of the elliptic
// problem.
//
// NOTES:
//
// OVERALL METHOD:
//
// 1. Input parameters (mass parameter) from stdin.
// 2. For each energy level H in the range, do
//    2.1 Input fixed points p=p3 and p4 in Delaunay coordinates, which are
//    on the section {g=0}.
//    2.2. Compute the complex integral $A^+$.
//    2.3. Compute the complex integral $B_in^{f,+}$.
// 3. For each energy level, we output one line to stdout:
//       H, Re(A^+(H)), Im(A^+(H)), Re(B_in^{f,+}), Im(B_in^{f,+}),
//       Re(B_in^{b,+}), Im(B_in^{b,+})

#include <stdio.h>
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE
#include <math.h>	// sqrt
#include <gsl/gsl_errno.h>      // gsl_set_error_handler_off
#include <rtbp.h> 	// DIM
#include "inner_ell.h" 	// re_inner_ell, im_inner_ell

int main( )
{
   double mu, H;

   double p3[DIM], p4[DIM];	// p=(l,L,g,G) fixed point in Delaunay

   double re_A, im_A;	// Re(A^+(H)), Im(A^+(H))

   double re_B_in_f, im_B_in_f;	// Re(A_in^{f,+}(H)), Im(A_in^{f,+}(H))
   double re_B_in_b, im_B_in_b;	// Re(A_in^{b,+}(H)), Im(A_in^{b,+}(H))

   // Input mass parameter from stdin.
   if(scanf("%le", &mu) < 1)
   {
      perror("main: error reading input");
      exit(EXIT_FAILURE);
   }

   // Stop GSL default error handler from aborting the program
   gsl_set_error_handler_off();

   // Input H, fixed point p
   while(scanf("%le %le %le %le %le %le %le %le %le", &H, 
	    p3, p3+1, p3+2, p3+3, p4, p4+1, p4+2, p4+3) == 9)
   {
      // Compute the complex integral $A^+$.
      // Integral goes from 0 to -14\pi.
      re_inner_ell(mu,14*M_PI,p3,&re_A);
      im_inner_ell(mu,14*M_PI,p3,&im_A);

      // Compute the complex integral $B_in^{f,+}$.
      // Integral goes from 0 to -12\pi.
      re_inner_ell(mu,12*M_PI,p4,&re_B_in_f);
      im_inner_ell(mu,12*M_PI,p4,&im_B_in_f);

      // Compute the complex integral $B_in^{b,+}$.
      // Integral goes from 0 to -2\pi.
      re_inner_ell(mu,2*M_PI,p3,&re_B_in_b);
      im_inner_ell(mu,2*M_PI,p3,&im_B_in_b);

      // For each energy level, we output one line to stdout:
      //    H, Re(A^+(H)), Im(A^+(H)), 
      //    Re(B_in^{f,+}), Im(B_in^{f,+}), 
      //    Re(B_in^{b,+}), Im(B_in^{b,+}), 
      if(printf("%e %.15e %.15e %.15e %.15e %.15e %.15e\n", H, re_A, im_A, 
	       re_B_in_f, im_B_in_f, 
	       re_B_in_b, im_B_in_b)<0)
      {
	 perror("main: error writting output");
	 exit(EXIT_FAILURE);
      }
   }
   exit(EXIT_SUCCESS);
}
