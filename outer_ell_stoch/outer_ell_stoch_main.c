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
//    - stability
//       unst/st flag
// 1. Input parameters from stdin:
// 
//    - mass parameter "mu"
//
// 2. Process input table, where each line corresponds to an energy level,
// in the following way:
//
//    2.1. Input data from stdin
//    H (l_p, L_p, g_p=0, G_p) (l_h, L_h, g_h=0/pi, G_h) omega_neg M
//
//    2.2. Compute integrals $B^+$ and $C^+$ using numerical integration
//
//    2.3. Output data to stdout
//    H \re(B^+) \im(B^+) \re(C^+) \im(C^+)
//
// NOTES
// =====
// Homoclinic point is either on the {g=0} section (unst_br1, st_br1) or
// {g=\pi} section (unst_br2, st_br2).
// 
// $M$ is the number of {g=ct} Poincare iterates that we needed in the shooting
// method to obtain the homoclinic point $z$ from $z_u$. This is found by the
// formula integration_time / (2\pi), because $\dot g \approx -1$. 
//
// We use upper integration limit N = M.
//
// Recall that omega_pos= -omega_neg.

// Headers
#include <stdio.h>
#include <stdlib.h>     // EXIT_SUCCESS, EXIT_FAILURE
#include <math.h>     // M_PI
#include <rtbp.h>       // DIM
#include <approxint.h>  // stability_t

#include <utils_module.h>       // dblcpy

// re_B_stoch, im_B_stoch, re_C_stoch, im_C_stoch
#include "outer_ell_stoch.h" 	

int main( )
{

   // "stability" flag specifies wheather we want to use the unstable branch
   // (=0) or stable branch (=1) of the manifold
   int stability;

   double mu, H;
   double t;		// integration time to reach z from z_u
   int M,N;
   double p[DIM];	// periodic point
   double zu[DIM];	// homoclinic point z_u
   double zs[DIM];	// homoclinic point z_s
   double omega_pos;	// $\omega_+^j$
   double omega_neg;	// $\omega_-^j$

   double re_B_res, im_B_res;	// $\re(B^+)$, $\im(B^+)$
   double re_C_res, im_C_res;	// $\re(C^+)$, $\im(C^+)$

   // auxiliary variables
   int status;
   stability_t st;

   // Input parameters from stdin: mass mu
   if(scanf("%le %d", &mu, &stability) < 2)
   {
      perror("main: error reading input");
      exit(EXIT_FAILURE);
   }

   st = (stability==0 ? UNSTABLE : STABLE);

   if (st == UNSTABLE)
   {
	   // Process input table, where each line corresponds to an energy level.
	   while((scanf("%le", &H) == 1) && 
		 (scanf("%le %le %le %le", p, p+1, p+2, p+3) == 4) &&
			 (scanf("%le %le %le %le", zu, zu+1, zu+2, zu+3) == 4) && 
			 (scanf("%le %le", &omega_neg, &t) == 2))
	   {
			   // Compute M
			   M = t/(2*M_PI);

			   // No need for zs since we are not computing re/im_C_stoch.
			   // Instead of fetching zs from intersecs_st_SECg_br1.res, 
			   // we set it to the symmetric point of zu. 
			   //dblcpy(zs,zu,DIM);
			   //zs[0] = -zu[0];
			   //zs[2] = -zu[2];

			  N = M;
			  omega_pos= -omega_neg;

			  // Compute integrals $B^+$ and $C^+$ using numerical integration
			  // Numerically, we observe that: re_C = re_B, im_C = -im_B.

			  re_B_stoch(mu, p, zu, omega_pos, &re_B_res, M, N);
			  im_B_stoch(mu, p, zu, omega_pos, &im_B_res, M, N);
			  
			  //re_C_stoch(mu, p, zs, omega_neg, &re_C_res, M, N);
			  //im_C_stoch(mu, p, zs, omega_neg, &im_C_res, M, N);

			  // Output data to stdout
			  //    H \re(B^+) \im(B^+) \re(C^+) \im(C^+)
			  if(printf("%e %.15e %.15e %.15e %.15e\n", 
				   H, re_B_res, im_B_res, re_B_res, -im_B_res)<0)
			  //if(printf("%e %.15e\n", H, im_B_res)<0)
			  {
				 perror("main: error writting output");
				 exit(EXIT_FAILURE);
			  }
			  fflush(NULL);
	   }
   }
   else if (st == STABLE)
   {
	   // Process input table, where each line corresponds to an energy level.
	   while((scanf("%le", &H) == 1) && 
		 (scanf("%le %le %le %le", p, p+1, p+2, p+3) == 4) &&
			 (scanf("%le %le %le %le", zs, zs+1, zs+2, zs+3) == 4) && 
			 (scanf("%le %le", &omega_pos, &t) == 2))
	   {
			   // Compute M
			   M = -t/(2*M_PI);	// Recall: t is negative for st branch

			   // No need for zu since we are not computing re/im_B_stoch.
			   // Instead of fetching zu from intersecs_unst_SECg_br1.res, 
			   // we set it to the symmetric point of zs. 
			   //dblcpy(zu,zs,DIM);
			   //zu[0] = -zs[0];
			   //zu[2] = -zs[2];

			  N = M;
			  omega_neg= -omega_pos;

			  // Compute integrals $B^+$ and $C^+$ using numerical integration
			  // Numerically, we observe that: re_C = re_B, im_C = -im_B.

			  //re_B_stoch(mu, p, zu, omega_pos, &re_B_res, M, N);
			  //im_B_stoch(mu, p, zu, omega_pos, &im_B_res, M, N);
			  
			  re_C_stoch(mu, p, zs, omega_neg, &re_C_res, M, N);
			  im_C_stoch(mu, p, zs, omega_neg, &im_C_res, M, N);

			  // Output data to stdout
			  //    H \re(B^+) \im(B^+) \re(C^+) \im(C^+)
			  if(printf("%e %.15e %.15e %.15e %.15e\n", 
				   H, re_C_res, -im_C_res, re_C_res, im_C_res)<0)
			  //if(printf("%e %.15e\n", H, im_B_res)<0)
			  {
				 perror("main: error writting output");
				 exit(EXIT_FAILURE);
			  }
			  fflush(NULL);
	   }
   }
   exit(EXIT_SUCCESS);
}
