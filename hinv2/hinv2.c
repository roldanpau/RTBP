// ===========================
// Invert Hamiltonian equation
// ===========================
// FILE:          $RCSfile: hinv2.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2012-06-06 10:52:33 $
//
// PURPOSE:
//
// NOTES:
//
// OVERALL METHOD:
//
// FUNCTIONS
// =========
//
// hinv
//    Invert Hamiltonian equation
//       H(x,y,px,py)=H0,
//    solving for the unknown py.

#include <stdio.h>	// fprintf
#include <math.h>	// sqrt, fabs
#include <rtbp.h>	// DIM
#include <hinv.h>	// ERR_CPLX_ROOTS

// name OF FUNCTION: hinv
// CREDIT: 
//
// PURPOSE
// =======
// Consider the Hamiltonian $H$ of the RTBP in rotating coordinates
// \[ H(x,y,px,py) = 1/2 (px^2+py^2) + ypx - xpy - mu_1/r1 - mu_2/r2. \]
// Let the value of the Hamiltonian H=H0 be given. Supose we know the value
// of three coordinates, x,y,px.
// This procedure inverts the Hamiltonian equation
//    H(x,y,px,py)=H0,
// solving for the unknown py.
// Since the Hamiltonian equation is a quadratic equation in p_y, there are
// two possible roots.
// WE USE THE FIRST BRANCH.
//
// PARAMETERS
// ==========
// mu
//    mass parameter for the RTBP
// H
//    energy value
// p
//    point, 4 coordinates: p=(x,y,p_x,p_y). 
//    On input, p[0:2] holds the known values of x, y and p_x.
//    On return of the this function, p[3] holds the value of p_y such that
//    H(x,y,p_x,p_y)=H.
// 
// RETURN VALUE
// ============
// Returns a non-zero error code to indicate an error and 0 to indicate
// success.
// If an integration error is encountered, the function returns a non-zero
// value, and p is unmodified:
// 
// ERR_CPLX_ROOTS
//    There are no real roots. Negative discriminant yields complex roots.

int hinv2(double mu, double H,double p[DIM])
{
   // larger mass on the left of origin, smaller mass on the right
   double mu1 = mu;	
   double mu2 = 1.0-mu;
   double x=p[0], y=p[1], px=p[2];
   double a=0.5;
   double b=-x;
   double r1=sqrt((x-mu2)*(x-mu2)+y*y);
   double r2=sqrt((x+mu1)*(x+mu1)+y*y);
   double c=0.5*px*px + y*px - mu1/r1 - mu2/r2 - H;
   double disc = b*b-4.0*a*c;	// discriminant

   // Check that discriminant is not negative
   if(disc<0)
   {
      fprintf(stderr, "hinv: no real roots\n");
      return(ERR_CPLX_ROOTS);
   }

   // For 1:7 resonance, WE USE THE FIRST BRANCH.
   // For 3:1 resonance, WE USE THE SECOND BRANCH.
   //p[3]=(-b+sqrt(disc))/(2.0*a);
   p[3]=(-b-sqrt(disc))/(2.0*a);
   return 0;
}
