// ==========================================
// Continuation of the shift A_in wrt energy H
// ==========================================
// FILE:          $RCSfile: inner_ell_stoch_main.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2011-07-01 09:20:28 $
//
// PURPOSE:
// Consider the family of 3:1 resonant periodic orbits in the RTBP.
// Define the integral (a part of $B_in^j$ defined in the paper)
//    A_in := i \int_0^{2\pi} \frac{\Delta H_ell^{1,+}\circ \lambda}
//       {-1+\mu\partial_G \Delta H_circ\circ \lambda} e^{i\tilde\lambda}
//       d\sigma.
// Compute the periodic orbit for a given (wide) range of energies, and
// compute how the integral $A_in$ changes when we change the energy H.
//
// The function $B_in^j = \frac{1-e^{i\omega}}{1-e^{i2\pi\nu}} A_in$ is
// computed in outer_ell_stoch/B_in.c
//
// NOTES:
//
// OVERALL METHOD:
//
// 1. Input parameters (mass parameter) from stdin.
// 2. For each energy level H in the range, do
//    2.1 Input fixed point p in Delaunay coordinates, which is on the section
//    {g=0}.
//    2.3. Compute the complex integral $A_in$.
// 3. For each energy level, we output one line to stdout:
//       H, Re(A_in), Im(A_in),

#include <stdio.h>
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE
#include <math.h>	// sqrt
#include <gsl/gsl_errno.h>      // gsl_set_error_handler_off
#include <rtbp.h> 	// DIM
#include "inner_ell_stoch.h" 	// re_inner_ell_stoch, im_inner_ell_stoch

int main( )
{
   double mu, H;

   double p[DIM];	// p=(l,L,g,G) fixed point in Delaunay

   double re_A_in, im_A_in;	// Re(A_in(H)), Im(A_in(H))

   // Input mass parameter from stdin.
   if(scanf("%le", &mu) < 1)
   {
      perror("main: error reading input");
      exit(EXIT_FAILURE);
   }

   // Stop GSL default error handler from aborting the program
   gsl_set_error_handler_off();

   // Input H, fixed point p
   while(scanf("%le %le %le %le %le", &H, p, p+1, p+2, p+3) == 5)
   {
      // Compute the complex integral $A_in$.
      // Integral goes from 0 to 2\pi.
      re_inner_ell_stoch(mu,2*M_PI,p,&re_A_in);
      im_inner_ell_stoch(mu,2*M_PI,p,&im_A_in);

      // For each energy level, we output one line to stdout:
      //    H, Re(A_in), Im(A_in), 
      if(printf("%e %.15e %.15e\n", H, re_A_in, im_A_in)<0)
      {
         perror("main: error writting output");
         exit(EXIT_FAILURE);
      }
	  fflush(NULL);
   }
   exit(EXIT_SUCCESS);
}
