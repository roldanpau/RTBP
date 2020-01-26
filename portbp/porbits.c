/*! \file porbits.c
    \brief Continuation of the periodic orbit wrt energy H.
    
    Consider the family of p:q almost-resonant periodic orbits in the RTBP.
    Compute the periodic orbit for a given (wide) range of energies when we
    change the energy H.

    \remark The initial condition for the periodic orbit is taken at 
    the perihelion.

    $Author: roldan $
    $Date: 2013-02-25 09:53:04 $
*/

#include <stdio.h>
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE
#include <string.h>	// strcmp

#include <gsl/gsl_errno.h>      // gsl_set_error_handler_off
#include <rtbp.h> 		// DIM
#include <initcond.h>
#include <prtbp_2d.h>
#include <hinv.h>
#include "portbpsym.h"

/**
  Main program.

  Input params (stdin): mass parameter, Poincare section, p/q resonance,
  number of cuts, energy range [H_lo, H_hi], energy increment H_delta.

  Output params (stdout): H, T(H), x(H)
*/
int main( )
{
   double mu, H, H_lo, H_hi, H_delta, x, py;

   section_t sec;	// Poincare section

   double resonance; 	// resonance=p/q (this will be 7/1 or 1/3)

   // number of cuts with section (for 1:3 resonance, k=2)
   int k;
   
   // periodic point, 2d cartesian coordinates: z=(x,px)
   double z[2]; 	

   // periodic point, 4d cartesian coordinates: z2=(x,y,px,py)
   double z2[DIM]; 	

   //double p[DIM]; 	// periodic point, Delaunay coordinates: p=(l,L,g,G)
   double T;		// period of periodic orbit

   // auxiliary variables
   int status;
   char section_str[10];        // holds input string "SEC1", "SEC2" etc

   // Input mass parameter, Poincare section, p/q resonance, number of cuts, energy range,
   // energy increment from stdin.
   if(scanf("%le %s %le %d %le %le %le", &mu, section_str, &resonance, &k, &H_lo, &H_hi,
	    &H_delta) < 7)
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

   for(H=H_lo; H<=H_hi; H+=H_delta)
   {
      fprintf(stderr, "H = %.15e\n", H);

      // Compute initial condition x_0 corresponding to p:q res. p.o. in 2BP.
      status = initcond(H, resonance, &x, &py);

      // Compute the fixed point x corresponding to (almost) res. p.o. in
      // RTBP.

      // Refine trajectory of the RTBP which is close to a periodic orbit,
      // until a true periodic orbit is obtained.
      status = portbpsym(mu,sec,H,k,&x);
      if(status)
      {
	 fprintf(stderr, \
	       "main: unable to find periodic orbit up to desired accuracy\n");
	 // this is not really an error
	 // exit(EXIT_FAILURE);
      }

      // Compute the period T of the p.o.
      z[0] = x;
      z[1] = 0;	// px
      if(prtbp_2d(mu,sec,H,k,z,&T))
      {
	 fprintf(stderr, "main: error computing period\n");
	 exit(EXIT_FAILURE);
      }

      // Transform periodic point z to Delaunay variables p.
      /*
      z[0] = x;
      z[1] = 0;	//px

      if(cardel_2d(mu,H,z,py,p))
      {
	 fprintf(stderr, "main: error obtaining Delaunay coordinates\n");
	 exit(EXIT_FAILURE);
      }
      */

      // Obtain coordinate p_y
      z2[0] = x;
      z2[1] = 0;	// y
      z2[2] = 0;	// px
      hinv(mu,sec,H,z2);

      // Output one line to stdout: H, T(H), x(H)
      // For the moment, we don't output p(H)...
      if(printf("%.15e %.15e %.15e %.15e %.15e %.15e\n", H, 2*T, z2[0], z2[1], z2[2], z2[3])<0)
      {
	 perror("main: error writting output");
	 exit(EXIT_FAILURE);
      }
      fflush(NULL);
   }
   exit(EXIT_SUCCESS);
}
