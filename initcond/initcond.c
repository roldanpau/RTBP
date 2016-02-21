/*! \file initcond.c
    \brief initial condition of p/q resonant p.o. in 2BP

    $Author: roldan $
    $Date: 2013-03-13 09:03:57 $
*/

#include <math.h>	// pow

/**
  initial condition of p/q resonant p.o. in 2BP

  Consider the Asteroid - Sun/Jupiter two body problem (2BP), i.e. mass of
  Jupiter is zero.
  Given an energy value "H", this function computes an initial condition of
  Asteroid (x,y,v_x,v_y) corresponding to r=p/q ressonant periodic orbit in
  the 2BP, where p=period of asteroid, q=period of Jupiter.
 
  We impose that initial condition is at the perihelion.
  This implies that position of Asteroid is A=(x,0) and velocity is
  v=(0,v_y), where x>0 and p_y>0. Therefore, we output as initial condition
  only the components (x,v_y).
 
  \param[in] H 	energy value where we will look for periodic orbit
  \param[in] r 	r=p/q resonance, where p=period of asteroid, q=period of Jupiter
  \param[out] x On exit, it contains the x component of the initial condition
  \param[out] vy On exit, it contains the v_y component of the initial condition

  \returns
  a non-zero error code to indicate an error and 0 to indicate
  success.
  */

//
// NOTES
// =====
// We use Vadim's computations to obtain the initial condition (see the
// scanned notes "kaloshin.pdf").

int initcond(double H, double r, double *x, double *vy)
{
   double L=pow(r,1.0/3.0);		// L^2 = semi-major axis
   double G=-H-1/(2*L*L);	// angular momentum
   double GL = G/L;
   double e=sqrt(1-GL*GL);	// eccentricity

   // Suppose i.c. at perihelion: pos=(x,0), vel=(0,v_y)

   *x=G*G/(1+e);
   *vy=(1+e)/G;
   return(0);
}
