// =====================
// Function $B^{f,+}(I)$
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
// \[ B^{f,+}(I) = B_{out}^{f,+}(I) + 
//    B_{in}^{f,+}(I) e^{i\mu\omega_{out}^f(I)} \]
//
// where 
// \[ B_{out}^{f,+} = \mu B^+ - \mu C^+ = i 2\mu \Im(B^+). \]
//
// Thus
// \[ \Re(B^f) = \Re(B_in^f)\cos(\mu\omega_out^f) - 
//    \Im(B_in^f)\sin(\mu\omega_out^f), \\
//
//    \Im(B^f) = 2\mu\Im(B^+) + \Im(B_in^f)\cos(\mu\omega_out^f) + 
//    \Re(B_in^f)\sin(\mu\omega_out^f), \]
//
// with
//    \[ \omega_out^f = -2\omega_+$. \]
//
//
// OVERALL METHOD
// ==============
//
// 2. Process input table, where each line corresponds to an energy level,
// in the following way:
//
//    2.1. Input data from stdin
//    H, \Im(B^+), \Re(B_{in}^f), \Im(B_{in}^f), \omega_pos
//
//    2.2. Compute complex function $B^{f,+}(I)$ (real and imaginary parts).
//
//    2.3. Output data to stdout
//    H \re(B^{f,+}) \im(B^{f,+})
//
// NOTES
// =====

// remember that, since Marcel changed signs in the last version of the paper,
// our computed value of \omega_+ is now the negative of Marcel's \omega_+.


#include <stdio.h>
#include <stdlib.h>     // EXIT_SUCCESS, EXIT_FAILURE
#include <math.h>	// sin, cos

const double mu=0.95387536e-3;

int main( )
{
//    H, \Im(B^+), \Re(B_{in}^f), \Im(B_{in}^f), \omega_pos
//
   double H;
   double imBpos, reBin, imBin, omega_pos;

   double reB, imB;	// \re(B^{f,+}) \im(B^{f,+})

   double omega_out;

   // Process input table, where each line corresponds to an energy level.
   while(scanf("%le %le %le %le %le",
            &H, &imBpos, &reBin, &imBin, &omega_pos) == 5)
   {
      omega_out = -2*omega_pos;

      // Compute complex function $B^{f,+}$
      reB = reBin*cos(mu*omega_out) - imBin*sin(mu*omega_out);
      imB = 2*mu*imBpos + imBin*cos(mu*omega_out) + reBin*sin(mu*omega_out);

      // Output data to stdout
      //    H \re(B^{f,+}) \im(B^{f,+})
      if(printf("%e %.15e %.15e\n", H, reB, imB)<0)
      {
         perror("main: error writting output");
         exit(EXIT_FAILURE);
      }
   }
   exit(EXIT_SUCCESS);
}

