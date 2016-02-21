// ===========================
// Invert Hamiltonian equation
// ===========================
// FILE:          $RCSfile: hinv2_main.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2011-03-10 14:32:22 $
//
// PURPOSE
// =======
// Consider the Hamiltonian $H$ of the RTBP in rotating coordinates
// \[ H(x,y,px,py) = 1/2 (px^2+py^2) + ypx - xpy - (1-mu)/r1 - mu/r2. \]
// Let the value of the Hamiltonian H=H0 be given. Supose we know the value
// of three coordinates, x,y,px.
// This program inverts the Hamiltonian equation
//    H(x,y,px,py)=H0,
// solving for the unknown py.
// Since the Hamiltonian equation is a quadratic equation in p_y, there are
// two possible roots. The caller must choose one root by providing us with
// an approximation to p_y.
//
// NOTES
// =====
//
// OVERALL METHOD:
//
// 1. Input parameters from stdin:
// 
// mu
//    mass parameter for the RTBP
// H
//    energy value
// p
//    point, 4 coordinates: p=(x,y,p_x,p_y). 
//    p[0:2] holds the known values of x, y and p_x.
//    p[3] holds the approximate value of p_y.
//
// 2. Invert Hamiltonian equation, solving for p_y.
// 3. Output solution p_y such that H(x,y,p_x,p_y)=H.

#include <stdio.h>
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE
#include <rtbp.h>	// DIM
#include "hinv2.h"	// hinv

int main( )
{
   double mu, H;
   double p[DIM];
   int status;

   // Input mass parameter, energy value, and point from stdin.
   if(scanf("%le %le %le %le %le %le", &mu, &H, p, p+1, p+2, p+3) < 6)
   {
      perror("main: error reading input");
      exit(EXIT_FAILURE);
   }

   // Invert Hamiltonian equation, solving for p_y.
   status=hinv2(mu,H,p);
   if(status)
   {
      fprintf(stderr, \
	    "main: error inverting Hamiltonian equation\n");
      exit(EXIT_FAILURE);
   }

   // Output solution p_y to stdout.
   if(printf("%.15le\n", p[3])<0)
   {
      perror("main: error writting output");
      exit(EXIT_FAILURE);
   }
   exit(EXIT_SUCCESS);
}
