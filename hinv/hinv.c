// ===========================
// Invert Hamiltonian equation
// ===========================
// FILE:          $RCSfile: hinv.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2012-12-19 10:49:23 $
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
#include <prtbp.h>	// section_t

const int ERR_CPLX_ROOTS=1;

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
//
// NOTE!!! It is better to solve 
//    C(x,y,x',y') = -(x'^2 + y'^2) + 2\Omega(x,y)
// or equivalently
//    H(x,y,x',y') = 1/2 (x'^2 + y'^2) - 1/2 (x^2+y^2) -F(x,y)
// (see Szebehelly, page 351)
// for |x'|, and determine its sign by hand (depending on which direction we
// cross the section). See Canalias and Masdemont: "Homoclinic and
// heteroclinic transfer trajectories..."
//
// NOTES
// =====
// Since the Hamiltonian equation is a quadratic equation in p_y, there are
// two possible roots.
// For the 3:1 resonance, one can use ONLY THE FIRST BRANCH.
// For the 1:7 resonance, one can use ONLY THE SECOND BRANCH.
//
// PARAMETERS
// ==========
// mu
//    mass parameter for the RTBP
// sec
//    type of Poincare section (sec = SEC1 or SEC2).
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
// If an error is encountered, the function returns a non-zero value, and p
// is unmodified:
// 
// ERR_CPLX_ROOTS
//    There are no real roots. Negative discriminant yields complex roots.

int hinv(double mu, section_t sec, double H,double p[DIM])
{
   // larger mass on the left of origin, smaller mass on the right
   double mu1 = mu;	
   double mu2 = 1.0-mu;
   double x=p[0], y=p[1], px=p[2];
   double r1sqr=(x-mu2)*(x-mu2)+y*y;
   double r2sqr=(x+mu1)*(x+mu1)+y*y;
   double r1=sqrt(r1sqr);
   double r2=sqrt(r2sqr);

   double vx=px+y;	// x component of the velocity
   double vy;		// y component of the velocity

   double F = mu1/r1+mu2/r2;
   double disc = -vx*vx+(x*x+y*y)+2*F+2*H;

   // Check that discriminant is not negative
   if(disc<0)
   {
      fprintf(stderr, "hinv: no real roots\n");
      return(ERR_CPLX_ROOTS);
   }

   if(sec==SEC1)
      vy = sqrt(disc);	// in section SEC1, we want vy>0
   else if(sec==SEC2)
      vy = -sqrt(disc);	// in section SEC2, we want vy<0

   p[3]=vy+x;
   return 0;
}
