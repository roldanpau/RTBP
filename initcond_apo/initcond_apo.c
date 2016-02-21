// =============================================
// initial condition of p/q resonant p.o. in 2BP
// =============================================
// FILE:          $RCSfile: initcond_apo.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2012-06-06 11:05:39 $
//
// FUNCTIONS
// =========
//
// initcond
// --------
// Consider the Asteroid - Sun/Jupiter two body problem (2BP).
// Given an energy value "h", this function computes an initial condition of
// Asteroid (x,y,p_x,p_y) corresponding to r=p/q ressonant periodic orbit in
// the 2BP.
//
// We impose that initial condition is at the apohelion. This implies that
// position of Asteroid is A=(-x,0) and velocity (and thus momentum) is
// v=(0,-v_y). Therefore, we output as initial condition only the components
// (-x,-v_y).

#include <math.h>	// pow

// name OF FUNCTION: initcond_apo
// CREDIT: 
//
// PURPOSE
// =======
// Consider the Asteroid - Sun/Jupiter two body problem (2BP), i.e. mass of
// Jupiter is zero.
// Given an energy value "h", this function computes an initial condition of
// Asteroid (x,y,v_x,v_y) corresponding to r=p/q ressonant periodic orbit in
// the 2BP, where p=period of asteroid, q=period of Jupiter.
//
// We impose that initial condition is at the apohelion. This implies that
// position of Asteroid is A=(-x,0) and velocity is v=(0,-v_y). Therefore, we
// output as initial condition only the components (-x,-v_y).
//
// PARAMETERS
// ==========
// 
// H
//    energy value where we will look for periodic orbit
// r
//    r=p/q resonance, where p=period of asteroid, q=period of Jupiter
// x
//    On exit, it contains the x component of the initial condition
// vy
//    On exit, it contains the v_y component of the initial condition
// 
// RETURN VALUE
// ============
// Returns a non-zero error code to indicate an error and 0 to indicate
// success.
//
// NOTES
// =====
// We use Vadim's computations to obtain the initial condition (see the
// scanned notes "kaloshin.pdf").
//
// Coordinates A=(-x,0), v=(0,-v_y) correspond to nonrotating coordinates.
// In rotating, we would have A=(x,0) instead.
//
// CALLS TO: 

int initcond_apo(double H, double r, double *x, double *vy)
{
   double L=pow(r,1.0/3.0);		// L^2 = semi-major axis
   double G=-H-1/(2*L*L);	// angular momentum
   double GL = G/L;
   double e=sqrt(1-GL*GL);	// eccentricity

   // Suppose i.c. at apohelion: pos=(-x,0), vel=(0,-v_y)

   *x=-G*G/(1-e);
   *vy=-(1-e)/G;
   return(0);
}
