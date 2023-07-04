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
// \[ B_{out}^j = -B^+ + C^+ = -i 2 \Im(B^+), \]
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
// NOTE: We noticed that in the paper there was a mistake in the formulas: mu
// was included twice, both in $\Delta H_{ell}^{1,+}$, and factored in front of
// the integrals in $B_{in}^j$ and $B_{out}^j$. Thus I decide to remove the
// extra mu factored in front of the integrals in this file.
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
#include <approxint.h>  // stability_t

int main( )
{
   // "stability" flag specifies wheather we want to use the unstable branch
   // (=0) or stable branch (=1) of the manifold
   int stability;

   double mu;

//    H, T, \Im(B^+), \Re(B_{in}), \Im(B_{in}), \omega_neg
   double H, T;
   double imBpos, reBin, imBin, omega;

   double nu;   // \mu\nu = T (where T = period mod 2\pi)

   double reB, imB;	// \re(B^j) \im(B^j)
   double reB_in, imB_in;	// \re(B^j_in) \im(B^j_in)
   double imB_out;	// \im(B^j_out)

   double omega_j;

   // auxiliary variables
   double a, b, c, d;
   stability_t st;

   // Input parameters from stdin: mass mu, stability flag
   if(scanf("%le %d", &mu, &stability) < 2)
   {
      perror("main: error reading input");
      exit(EXIT_FAILURE);
   }

   st = (stability==0 ? UNSTABLE : STABLE);

   // Process input table, where each line corresponds to an energy level.
   while(scanf("%le %le %le %le %le %le",
            &H, &T, &imBpos, &reBin, &imBin, &omega) == 6)
   {
	   if (st == UNSTABLE)	// omega = omega_neg
	   {
		  omega_j = -2*omega;
	   }
	   else if (st == STABLE)	// omega = omega_pos
	   {
		  omega_j = 2*omega;
	   }

      nu = fmod(T, 2*M_PI) / mu;

      a = 1-cos(omega_j);
      b = sin(omega_j);
      c = 1-cos(2*M_PI*nu);
      d = sin(2*M_PI*nu);

      // Compute complex function $B^j_in$
      reB_in = 1/(c*c+d*d)*((a*c+b*d)*reBin - (a*d-b*c)*imBin);
      imB_in = 1/(c*c+d*d)*((a*d-b*c)*reBin + (a*c+b*d)*imBin);

      // Compute complex function $B^j_out$
      imB_out = -2*imBpos;

      // Compute complex function $B^j$
      reB = reB_in;
      imB = imB_out + imB_in;

      // Output data to stdout
      //    H \re(B^j) \im(B^j)
      if(printf("%e %.15e %.15e %.15e %.15e %.15e\n", H, reB_in, imB_in,
                  imB_out, reB, imB)<0)
      {
         perror("main: error writting output");
         exit(EXIT_FAILURE);
      }
   }
   exit(EXIT_SUCCESS);
}

