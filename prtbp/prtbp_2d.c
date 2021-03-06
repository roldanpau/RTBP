/*! \file
    \brief Poincare map of the Restricted Three Body Problem (2D)

    $Author: roldan $
    $Date: 2013-03-26 22:22:14 $
*/

#include <stdio.h>	// fprintf
#include <stdlib.h>	// EXIT_FAILURE
#include <math.h>	// sqrt, fabs
#include <rtbp.h>	// DIM, rtbp
#include <hinv.h>	// hinv

#include <section.h>
#include <prtbp_nl_2d_module.h>
#include "prtbp.h"	// prtbp, prtbp_inv

int prtbp_2d(double mu, section_t sec, double H, int cuts, double p[2], double *ti)
{
	return(prtbp_nl_2d(mu,sec,H,cuts,p,ti));

   double x[DIM];

   // auxiliary variables
   int result;

   x[0]=p[0];	// x
   x[1]=0.0;	// y
   x[2]=p[1]; 	// p_x

   // Compute x[3]=p_y by inverting the Hamiltonian.
   result = hinv(mu,sec,H,x);
   if(result)
   {
      fprintf(stderr, "prtbp_2d: error inverting the Hamiltonian\n");
      return(1);
   }

   // Compute Poincare map (and integration time).
   // On exit, "prtbp" guarantees that x is precisely on the Poincare section.
   if(prtbp(mu,sec,cuts,x,ti))
   {
      fprintf(stderr, "prtbp_2d: error computing poincare map\n");
      return(1);
   }
   // Set the image point $P^n(p)$.
   p[0]=x[0]; 	// x'
   p[1]=x[2];	// p_x'
   return(0);
}

int prtbp_2d_inv(double mu, section_t sec, double H, int cuts, double p[2], double *ti)
{
   double x[DIM];

   // auxiliary variables
   int result;

   x[0]=p[0];	// x
   x[1]=0.0;	// y
   x[2]=p[1]; 	// p_x

   // Compute x[3]=p_y by inverting the Hamiltonian.
   result = hinv(mu,sec,H,x);
   if(result)
   {
      fprintf(stderr, "prtbp_2d_inv: error inverting the Hamiltonian\n");
      return(1);
   }

   // Compute inverse Poincare map (and integration time).
   // On exit, "prtbp_inv" guarantees that x is precisely on the Poincare section.
   if(prtbp_inv(mu,sec,cuts,x,ti))
   {
      fprintf(stderr, "prtbp_2d_inv: error computing inverse poincare map\n");
      return(1);
   }
   // Set the image point $P^{-n}(p)$.
   p[0]=x[0]; 	// x'
   p[1]=x[2];	// p_x'
   return(0);
}
