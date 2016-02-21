// =====================
// Function $B^{b,+}(I)$
// =====================
// FILE:          $RCSfile: B_b.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2011-07-22 08:48:27 $
//
// PURPOSE
// =======
// This program computes the complex function
// 
// \[ B^{b,+}(I) = B_{in}^{b,+}(I) + 
//    B_{out}^{b,+}(I) e^{i\mu\omega_{in}^b(I)} \]
//
// with
// 
// \[ \Re(B^b) = \Re(B_in^b) - 2\mu\Im(B^+) \sin(\mu\omega_in^b), \\
// \[ \Im(B^b) = \Im(B_in^b) + 2\mu\Im(B^+) \cos(\mu\omega_in^b). \\
//
//
// OVERALL METHOD
// ==============
//
// 2. Process input table, where each line corresponds to an energy level,
// in the following way:
//
//    2.1. Input data from stdin
//    H, \Im(B^+), \Re(B_{in}^b), \Im(B_{in}^b), \omega_{in}^b
//
//    2.2. Compute complex function $B^{b,+}(I)$ (real and imaginary parts).
//
//    2.3. Output data to stdout
//    H \re(B^{b,+}) \im(B^{b,+})
//
// NOTES
// =====


#include <stdio.h>
#include <stdlib.h>     // EXIT_SUCCESS, EXIT_FAILURE
#include <math.h>	// sin, cos

const double mu=0.95387536e-3;

int main( )
{
//    H, \Im(B^+), \Re(B_{in}^b), \Im(B_{in}^b), \omega_{in}^b
   double H;
   double imBpos, reBin, imBin, omega_in;

   double reB, imB;	// \re(B^{b,+}) \im(B^{b,+})

   // Process input table, where each line corresponds to an energy level.
   while(scanf("%le %le %le %le %le",
            &H, &imBpos, &reBin, &imBin, &omega_in) == 5)
   {
      // Compute complex function $B^{b,+}$
      reB = reBin - 2*mu*imBpos*sin(mu*omega_in);
      imB = imBin + 2*mu*imBpos*cos(mu*omega_in);

      // Output data to stdout
      //    H \re(B^{b,+}) \im(B^{b,+})
      if(printf("%e %.15e %.15e\n", H, reB, imB)<0)
      {
         perror("main: error writting output");
         exit(EXIT_FAILURE);
      }
   }
   exit(EXIT_SUCCESS);
}

