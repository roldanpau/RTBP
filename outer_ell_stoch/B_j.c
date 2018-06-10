// =====================
// Function $B^j(I)$
// =====================
// FILE:          $RCSfile: B_f.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2011-07-22 08:48:27 $
//
// PURPOSE
// =======
// This program computes the complex function
// 
// \[ B^j(I) = 
//       B_{in}^j(I) + B_{out}^j(I) \]
//
// where 
// \[ B_{in}^j = \frac{1-e^{i\omega^j}}{1-e^{i2\pi\nu}} A_in, \]
// \[ B_{out}^j = -\mu B^+ + \mu C^+ = -i 2\mu \Im(B^+), \]
// \[ \nu = T/mu. \]
//
// Thus
// \[ \Re(B^j) = ..., \\
//    \Im(B^j) = ..., \]
// (see my notes)
//
// with
//    \[ \omega^j = -2\omega_-^j. \]
//
//
// OVERALL METHOD
// ==============
//
// 2. Process input table, where each line corresponds to an energy level,
// in the following way:
//
//    2.1. Input data from stdin
//    H, T, \Im(B^+), \Re(B_{in}), \Im(B_{in}), \omega_neg
//
//    2.2. Compute complex function $B^j(I)$ (real and imaginary parts).
//
//    2.3. Output data to stdout
//    H \re(B^j) \im(B^j)


#include <stdio.h>
#include <stdlib.h>     // EXIT_SUCCESS, EXIT_FAILURE
#include <math.h>	// sin, cos

const double mu=0.95387536e-3;

int main( )
{
//    H, T, \Im(B^+), \Re(B_{in}), \Im(B_{in}), \omega_neg
   double H, T;
   double imBpos, reBin, imBin, omega_neg;

   double nu;   // \mu\nu = T (where T = period mod 2\pi)

   double reB, imB;	// \re(B^j) \im(B^j)

   double omega_j;

   // auxiliary variables
   double a, b, c, d;

   // Process input table, where each line corresponds to an energy level.
   while(scanf("%le %le %le %le %le %le",
            &H, &T, &imBpos, &reBin, &imBin, &omega_neg) == 6)
   {
      omega_j = -2*omega_neg;
      nu = fmod(T, 2*M_PI) / mu;

      a = 1-cos(omega_j);
      b = sin(omega_j);
      c = 1-cos(2*M_PI*nu);
      d = sin(2*M_PI*nu);

      // Compute complex function $B^j$
      reB = 1/(c*c+d*d)*((a*c+b*d)*reBin - (a*d-b*c)*imBin);
      imB = -2*mu*imBpos + 1/(c*c+d*d)*((a*d-b*c)*reBin + (a*c+b*d)*imBin);

      // Output data to stdout
      //    H \re(B^j) \im(B^j)
      if(printf("%e %.15e %.15e\n", H, reB, imB)<0)
      {
         perror("main: error writting output");
         exit(EXIT_FAILURE);
      }
   }
   exit(EXIT_SUCCESS);
}

