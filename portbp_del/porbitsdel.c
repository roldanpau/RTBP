// ===============================================
// Continuation of the periodic orbit wrt energy H
// ===============================================
// FILE:          $RCSfile: porbitsdel.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2012-06-06 12:59:19 $
//
// PURPOSE:
// Consider the family of p:q almost-resonant periodic orbits in the RTBP.
// Compute the periodic orbit for a given (wide) range of eccentricities when
// we change the eccentricity e.
//
// NOTES:
// We Look for a periodic orbit as a fixed point of the Poincare map. 
// As an initial condition to this fixed point, we can choose a point at the
// perihelion or the apohelion (both are on the symmetry line). 
// We choose initial cond. at the PERIHELION.
//
// OVERALL METHOD:
//
// 1. Input parameters (mass parameter, Poincare section, number of cuts,
//    eccentricity range, eccentricity increment) from stdin.
// 2. For each eccentricity e in the range, do
//    2.1. Compute initial condition z_0 = (g_0=0, G_0) corresponding to p:q
//         res. p.o. in 2BP.
//    2.2. Compute the fixed point z = (g, G) corresponding to almost-res.
//         p.o. in RTBP.
//    2.3. Compute the period T of the p.o.
//    2.4. Invert hamiltonian to get the L component of fixed point.
//    2.4. Output one line to stdout:
//            e, H(e), T(e), z(e)=(l,L,g,G)

#include <stdio.h>
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE
#include <math.h>	// pow, M_PI

#include <gsl/gsl_errno.h>      // gsl_set_error_handler_off
#include <rtbp.h> 		// DIM
#include <prtbpdel.h>		// section_t
#include <prtbpdel_2d.h>
#include "portbpdel.h"

int main( )
{
   double L0 = pow(3,-1.0/3.0);	// for 3:1 resonance

   double mu, e, e_lo, e_hi, e_delta, g, G;

   section_t sec;	// Poincare section

   double H;		// energy of initial condition
   int k;		// number of cuts with section
   int status;

   double z[2]; 	// periodic point, Delaunay coordinates: (g,G)
   double p[DIM];	// periodic point, Delaunay coordinates: (l,L,g,G) 
   double T;		// period of periodic orbit
   double L;		// L component of fixed point

   // auxiliary variables
   char section_str[10];        // holds input string "SEC1", "SEC2" etc
   double z_aux[2]; 		// periodic point, Delaunay coordinates: (g,G)

   // Input mass parameter, Poincare section, number of cuts, energy range,
   // energy increment from stdin.
   if(scanf("%le %s %d %le %le %le", &mu, section_str, &k, &e_lo, &e_hi, &e_delta) < 6)
   {
      perror("main: error reading input");
      exit(EXIT_FAILURE);
   }

   if (strcmp(section_str,"SEC1") == 0)
      sec = SEC1;
   else if (strcmp(section_str,"SEC2") == 0)
      sec = SEC2;
   else
   {
      perror("main: error reading section string");
      exit(EXIT_FAILURE);
   }

   // Stop GSL default error handler from aborting the program
   gsl_set_error_handler_off();

   for(e=e_lo; e<=e_hi; e+=e_delta)
   {

      // Compute initial condition x_0 corresponding to p:q res. p.o. in 2BP.
      g = 0.0;
      G = L0*sqrt(1.0-e*e);
      H = -1.0/(2*L0*L0)-G;

      // Compute the fixed point x corresponding to almost-res. p.o. in RTBP.
      z[0] = g;
      z[1] = G;

      // Refine trajectory of the RTBP which is close to a periodic orbit,
      // until a true periodic orbit is obtained.
      status = portbp_del(mu,sec,H,k,z);
      if(status)
      {
	 fprintf(stderr, \
	       "main: unable to find periodic orbit up to desired accuracy\n");
	 // this is not really an error
	 // exit(EXIT_FAILURE);
      }
      // Periodic orbit is actually on symmetry axis, so we impose g=0. 
      z[0]=0;

      // Compute the period T of the p.o.
      z_aux[0]=z[0];
      z_aux[1]=z[1];
      if(prtbp_del_2d(mu,sec,H,k,z_aux,&T))
      {
	 fprintf(stderr, "main: error computing period\n");
	 exit(EXIT_FAILURE);
      }

      // Invert hamiltonian to get the L component of fixed point
      switch(sec)
      {
	 case SEC1 :
	    {
	       p[0] = 0.0;         // l=0
	       break;
	    }
	 case SEC2 :
	    {
	       p[0] = M_PI;        // l=pi
	       break;
	    }
      }
      p[1]=L0;		// approximate L
      p[2]=z[0]; 	// g
      p[3]=z[1];	// G
      if(hinv_del(mu,H,p))
      {
	 fprintf(stderr, "main: error inverting hamiltonian\n");
	 exit(EXIT_FAILURE);
      }

      // Output one line to stdout: e, H(e), T(e), p(e)
      if(printf("%e %.15e %.15e %.15e %.15e %.15e %.15e\n", e, H, T, 
	       p[0], p[1], p[2], p[3])<0)
      {
	 perror("main: error writting output");
	 exit(EXIT_FAILURE);
      }
   }
   exit(EXIT_SUCCESS);
}
