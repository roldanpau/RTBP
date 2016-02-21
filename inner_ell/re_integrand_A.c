// ==========================================
// Continuation of the shift A^+ wrt energy H
// ==========================================
// FILE:          $RCSfile: re_integrand_A.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2011-03-24 14:06:08 $
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
// 1. Input parameters (mass parameter, energy range, energy increment) from
//    stdin.
// 2. For each energy level H in the range, do
//    2.1 Compute initial condition x_0 corresponding to 7:1 res. p.o. in
//        2BP.
//    2.2 Compute the fixed point x corresponding to almost-res. p.o. in
//        RTBP.
//    2.3 Transform fixed point from euclidean to Delaunay coordinates.
//    2.4. Compute the complex integral $A^+$.
// 3. For each energy level, we output one line to stdout:
//       H, x(H), Re(A^+(H)), Im(A^+(H)), modulus(A^+(H))

#include <stdio.h>
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE
#include <math.h>	// sqrt
#include <gsl/gsl_errno.h>      // gsl_set_error_handler_off
#include <initcond.h>
#include <portbpsym.h>
#include <rtbp.h> 	// DIM
#include "inner_ell.h" 	// re_inner_ell, im_inner_ell

struct iparams_inner_ell
{
   double mu;
   double x[DIM];
};

int main( )
{
   double mu, H, H_lo, H_hi, H_delta, x, py;

   double z[2];		// z=(x,px) fixed point
   double Y[DIM];	// Y=(l,L,g,G) fixed point in Delaunay
   int k=6;		// number of cuts with section {y=0}

   double re_A, im_A, mod_A;	// Re(A^+(H)), Im(A^+(H)), modulus(A^+(H))
   int status;

   double s, re_A_int;

   // Input mass parameter, energy range, energy increment from stdin.
   if(scanf("%le %le %le %le", &mu, &H_lo, &H_hi, &H_delta) < 4)
   {
      perror("main: error reading input");
      exit(EXIT_FAILURE);
   }

   // Stop GSL default error handler from aborting the program
   gsl_set_error_handler_off();

   H=H_lo;

//    2.1 Compute initial condition x_0 corresponding to 7:1 res. p.o. in
//    2BP.
   status = initcond(H, 7, &x, &py);

//    2.2 Compute the fixed point x corresponding to almost-res. p.o. in
//    RTBP.

   // Refine trajectory of the RTBP which is close to a periodic orbit, until
   // a true periodic orbit is obtained.
   status = portbpsym(mu,H,k,&x,py);
   if(status)
   {
      fprintf(stderr, \
	    "main: unable to find periodic orbit up to desired accuracy\n");
      // this is not really an error
      // exit(EXIT_FAILURE);
   }

//    2.3 Transform fixed point from euclidean to Delaunay coordinates.
   z[0] = x;
   z[1] = 0;	// px
   status = cardel_2d(mu,H,z,py,Y);
   if(status)
   {
      fprintf(stderr, "main: error obtaining Delaunay coordinates\n");
      exit(EXIT_FAILURE);
   }

   struct iparams_inner_ell params;
   params.mu = mu;
   params.x[0] = Y[0];
   params.x[1] = Y[1];
   params.x[2] = Y[2];
   params.x[3] = Y[3];

   for(s=0; s>-14*M_PI; s-= (M_PI/10.0))
   {
      re_A_int = re_integrand_inner_ell(s,&params);
      printf("%e %.15e\n", s, re_A_int);
   }
   exit(EXIT_SUCCESS);
}
