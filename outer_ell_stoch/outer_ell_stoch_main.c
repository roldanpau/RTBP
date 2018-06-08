// =================================
// Outer Map of the Elliptic Problem
// =================================
// FILE:          $RCSfile: outer_ell_main.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2011-07-22 08:47:46 $
//
// PURPOSE
// =======
// This program computes the shift of the outer map of the elliptic problem,
// given by the integral
// 
// \[ \Omega^+(I) = -\mu B^+(I) +\mu C^+(I), \]
//
// where 
// 
// \[ B^+(I) = \int_0^{14N\pi} 
//    f(\gamma_h(s)) e^{it(s)} - f(\gamma_p(s)) e^{i(t(s)-\omega_0^+)} ds, \]
// \[ C^+(I) = \int_0^{-14N\pi}
//    f(\gamma_h(s)) e^{it(s)} - f(\gamma_p(s)) e^{i(t(s)+\omega_0^-)} ds. \]
//
// Here, $\gamma_h(s)$ is the homoclinic trajectory, and $\gamma_p(s)$ is the
// periodic trajectory.
// (see notes "Inner and outer dynamics" by Marcel).
//
// OVERALL METHOD
// ==============
//
// 1. Input parameters from stdin:
// 
//    - mass parameter "mu"
//
// 2. Process input table, where each line corresponds to an energy level,
// in the following way:
//
//    2.1. Input data from stdin
//    H (l_p, L_p, g_p=0, G_p) (l_h, L_h, g_h=0, G_h) omega_neg M
//
//    2.2. Compute integrals $B^+$ and $C^+$ using numerical integration
//
//    2.3. Output data to stdout
//    H \re(B^+) \im(B^+) \re(C^+) \im(C^+)
//
// NOTES
// =====
// $M$ is the number of iterates that we needed in the shooting method to
// obtain the homoclinic point (program approxint).
//
// We use upper integration limit N = M.
//
// Recall that omega_pos= -omega_neg.

// Headers
#include <stdio.h>
#include <stdlib.h>     // EXIT_SUCCESS, EXIT_FAILURE
#include <math.h>     // M_PI
#include <rtbp.h>       // DIM

#include <utils_module.h>       // dblcpy

// re_B_stoch, im_B_stoch, re_C_stoch, im_C_stoch
#include "outer_ell_stoch.h" 	

int main( )
{
   double mu, H;
   int M,N;
   double p[DIM];	// periodic point
   double zu[DIM];	// homoclinic point z_u
   double zs[DIM];	// homoclinic point z_s
   double omega_pos;	// $\omega_+^j$
   double omega_neg;	// $\omega_-^j$

   double re_B_res, im_B_res;	// $\re(B^+)$, $\im(B^+)$
   double re_C_res, im_C_res;	// $\re(C^+)$, $\im(C^+)$

   // auxiliary variables
   double t;
   int status;

   // Input parameters from stdin: mass mu
   if(scanf("%le", &mu) < 1)
   {
      perror("main: error reading input");
      exit(EXIT_FAILURE);
   }

   // Process input table, where each line corresponds to an energy level.
   /*
   while((scanf("%le", &H) == 1) && 
	 (scanf("%le %le %le %le", p, p+1, p+2, p+3) == 4) &&
         (scanf("%le %le %le %le", zu, zu+1, zu+2, zu+3) == 4) && 
         (scanf("%le %le %le %le", zs, zs+1, zs+2, zs+3) == 4) && 
         (scanf("%le %d", &omega_neg, &M) == 2))
	 */
   while((scanf("%le", &H) == 1) && 
	 (scanf("%le %le %le %le", p, p+1, p+2, p+3) == 4) &&
         (scanf("%le %le %le %le", zu, zu+1, zu+2, zu+3) == 4) && 
         (scanf("%le %d", &omega_neg, &M) == 2))
   {
       // Instead of fetching zs from intersecs_st_SECg_br1.res, 
       // we set it to the symmetric point of zu. 
       dblcpy(zs,zu,DIM);
       zs[0] = -zu[0];
       zs[2] = -zu[2];

      N = M;
      omega_pos= -omega_neg;

      // Compute integrals $B^+$ and $C^+$ using numerical integration
      // Numerically, we observe that: re_C = re_B, im_C = -im_B.

      re_B_stoch(mu, p, zu, omega_pos, &re_B_res, M, N);
      im_B_stoch(mu, p, zu, omega_pos, &im_B_res, M, N);
      
      re_C_stoch(mu, p, zs, omega_neg, &re_C_res, M, N);
      im_C_stoch(mu, p, zs, omega_neg, &im_C_res, M, N);

      // Output data to stdout
      //    H \re(B^+) \im(B^+) \re(C^+) \im(C^+)
      if(printf("%e %.15e %.15e %.15e %.15e\n", 
	       H, re_B_res, im_B_res, re_C_res, im_C_res)<0)
      //if(printf("%e %.15e\n", H, im_B_res)<0)
      {
         perror("main: error writting output");
         exit(EXIT_FAILURE);
      }
      fflush(NULL);
   }
   exit(EXIT_SUCCESS);
}
