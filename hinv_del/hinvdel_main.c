// ===========================
// Invert Hamiltonian equation
// ===========================
// FILE:          $RCSfile: hinvdel_main.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2012-05-17 12:41:46 $
//
// PURPOSE
// =======
// Consider the Hamiltonian $H$ of the RTBP in rotating Delaunay coordinates
// \[ H(l,L,g,G) = -1/(2L^2) -G + \mu \Delta H_{circ}. \]
// Let the value of the Hamiltonian H=H0 be given. Supose we know the value
// of three coordinates, l,g,G.
// This procedure inverts the Hamiltonian equation
//    H(l,L,g,G)=H0,
// solving for the unknown L.
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
//    point, 4 coordinates: p=(l,L,g,G).
//    p[0],p[2],p[3] holds the known values of l, g and G.
//
// 2. Invert Hamiltonian equation, solving for L.
// 3. Output solution L such that H(l,L,g,G)=H.

#include <stdio.h>
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE
#include <rtbp.h>	// DIM
#include "hinvdel.h"	// hinv_del

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

   // Invert Hamiltonian equation, solving for L.
   status=hinv_del(mu,H,p);
   if(status)
   {
      fprintf(stderr, \
	    "main: error inverting Hamiltonian equation\n");
      exit(EXIT_FAILURE);
   }

   // Output solution L to stdout.
   if(printf("%.15le\n", p[1])<0)
   {
      perror("main: error writting output");
      exit(EXIT_FAILURE);
   }
   exit(EXIT_SUCCESS);
}
