/*! \file cardel.h
    \brief Cartesian to Delaunay Change of Coordinates

    Obtain rotating Delaunay coordinates (l,L,g,G) from rotating cartesian
    (x,y,px,py).

    $Author: roldan $
    $Date: 2013-03-26 22:15:14 $
*/

#include <rtbp.h>	//DIM

/**
  Obtain rotating Delaunay coordinates (l,L,g,G) from rotating cartesian
  (x,y,px,py).

  We have to choose the sign of $u$ (eccentric anomaly) correctly.
  For this, let us consider the (forward) trajectory through the point
  $(x,y,p_x,p_y)$.
  We consider \f$\dot r\f$ (which is known). 
  Since \f$\dot \ell>0\f$, Asteroid moves in a counter-clockwise fashion.
  Then, if \f$\dot r>0\f$ means that Asteroid is leaving the perihelion, so
  $u>0$. Conversely, if \f$\dot r<0\f$ means that Asteroid is approaching the
  perihelion, so $u<0$.
 
  On output, we normalize the angles $l,g$ between \f$[0,2\pi)\f$.
 
  Delaunay variables are 2BP variables, so they assume \f$\mu=0\f$.
  Thus this function does not change if we are using the "big primary to the
  left" or the "big primary to the right" convention.

  \param[in] X point in rotating cartesian coordinates, X=(x,y,p_x,p_y).
  \param[out] Y point in rotating Delaunay coordinates, Y=(l,L,g,G).
  \return void
*/    
void cardel(double X[DIM], double Y[DIM]);
