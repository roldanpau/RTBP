// =================================
// Outer Map of the Elliptic Problem
// =================================
// FILE:          $RCSfile: re_integrand_B.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2011-06-13 19:50:20 $
//
// PURPOSE
// =======
// This program computes the shift of the outer map of the elliptic problem,
// given by the integral
// 
// \[ \Omega^+(I) = -B^+(I) - C^+(I), \]
//
// where 
// 
// \[ B^+(I) = -\mu \int_0^{14M\pi} 
//    f(\gamma_h(s)) e^{it(s)} - f(\gamma_p(s)) e^{i(t(s)+\omega_0^+)} ds, \]
// \[ C^+(I) = -\mu \int_{-14M\pi}^0
//    f(\gamma_h(s)) e^{it(s)} - f(\gamma_p(s)) e^{i(t(s)-\omega_0^-)} ds. \]
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
//    - integration limit "M"
//
// 2. Process input table, where each line corresponds to an energy level,
// in the following way:
//
//    2.1. Input data from stdin
//    H (l_p, L_p, g_p=0, G_p) (l_h, L_h, g_h=0, G_h)^+ 
//    (l_h, L_h, g_h=0, G_h)^- omega_pos omega_neg
//
//    2.2. Compute integrals $B^+$ and $C^+$ using numerical integration
//
//    2.3. Output data to stdout
//    H \re(B^+) \im(B^+) \re(C^+) \im(C^+)
//
// NOTES
// =====
// Integration limit $M$ is fixed for all level of energies.
// In theory, we could choose an optimal $M$ for each level of energy by
// looking at the results of program omega_test (see my notes splitting.pdf).

// Headers
#include <stdio.h>
#include <stdlib.h>     // EXIT_SUCCESS, EXIT_FAILURE
#include <math.h>	// M_PI
#include <rtbp.h>       // DIM
#include "outer_ell.h" 	// re_B_pos, im_B_pos, re_C_pos, im_C_pos

struct iparams_outer_ell
{
   double mu;
   double p[DIM];
   double h[DIM];
   double omega;	// either $\omega_0^+$ or $\omega_0^-$
};

struct iparams_outer_ell_1
{
   double mu;
   double h[DIM];
};

struct iparams_outer_ell_2
{
   double mu;
   double p[DIM];
   double omega;	// either $\omega_0^+$ or $\omega_0^-$
};

int main( )
{
   double mu, H;
   int M;
   double p[DIM];	// periodic point
   double h_pos[DIM];	// homoclinic point h^+
   double h_neg[DIM];	// homoclinic point h^-
   double omega_pos;	// $\omega_0^+$
   double omega_neg;	// $\omega_0^-$

   double s;
   double re_B_int, re_B1_int, re_B2_int;

   // Input parameters from stdin
   if(scanf("%le %d", &mu, &M) < 2)
   {
      perror("main: error reading input");
      exit(EXIT_FAILURE);
   }

   // Process input table, where each line corresponds to an energy level.
   scanf("%le %le %le %le %le %le %le %le %le %le %le %le %le %le %le", 
	    &H, 
	    p, p+1, p+2, p+3, 
	    h_pos, h_pos+1, h_pos+2, h_pos+3, 
	    h_neg, h_neg+1, h_neg+2, h_neg+3, 
	    &omega_pos, &omega_neg);

   struct iparams_outer_ell params;
   struct iparams_outer_ell_1 params1;
   struct iparams_outer_ell_2 params2;

   params.mu = mu;
   params.h[0] = h_pos[0];
   params.h[1] = h_pos[1];
   params.h[2] = h_pos[2];
   params.h[3] = h_pos[3];
   params.p[0] = p[0];
   params.p[1] = p[1];
   params.p[2] = p[2];
   params.p[3] = p[3];
   params.omega = omega_pos;

   params1.mu = mu;
   params1.h[0] = h_pos[0];
   params1.h[1] = h_pos[1];
   params1.h[2] = h_pos[2];
   params1.h[3] = h_pos[3];

   params2.mu = mu;
   params2.p[0] = p[0];
   params2.p[1] = p[1];
   params2.p[2] = p[2];
   params2.p[3] = p[3];
   params2.omega = omega_pos;

   for(s=0; s<= 14*M*M_PI; s+= (M_PI/4.0))
   {
      re_B1_int = re_integrand_B1(s,&params1);
      re_B2_int = re_integrand_B2(s,&params2);
      re_B_int = re_integrand_B(s,&params);
      printf("%e %.15e %.15e %.15e\n", s, re_B1_int, re_B2_int, re_B_int);
   }
   exit(EXIT_SUCCESS);
}
