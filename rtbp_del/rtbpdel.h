/*! \file
    \brief Restricted Three Body Problem equations in Delaunay coordinates

    $Author: roldan $
    $Date: 2013-03-26 22:26:03 $
*/

#define DIM 4	///< dimension of the (planar) RTBP
#define ERR_COLLISION 1
double Hamilt_del(double mu, const double *x);

// Compute the eccentric anomaly u.
// We need to solve u-e sin(u) = l for u in (0,2*pi).
double eccentric(double e, double l);

/**
  Restricted Three Body Problem equations in Delaunay coordinates

  \authors Marcel Guardia and Pau Roldan

  The function computes the vectorfield of the RTBP at a point named "x".
  Specifically, we use the equations of motion of the planar circular
  RTBP in Delaunay coordinates (l,L,g,G), in the rotating (synodic)
  coordinate system.

  We follow the convention to place the large mass (Sun) to the left of
  the origin, and the small mass (Jupiter) to the right.
  This is opposite to the usual astrodynamics convention.

  This system is Hamiltonian, with Hamiltonian function
  \f[   H(l,L,g,G) = -1/(2L^2) - G + R, \f]
  where R is the perturbation relative to the two body problem.
  See Szebehelly, page 364.
  See Marcel's notes "Inner and outer dynamics".
  See also my hand-written notes "RTBP vector field in Delaunay".
  The function sets the value of the vectorfield in the variable "y".
  It returns a status code (success/error).

  We normalize $\ell$ between 0 and 2\pi to compute the vectorfield.
  This is done so that function "eccentric" is more precise.

  \param[in] t
  adimensional time at which the vectorfield is evaluated. Since this is
  an autonomous ODE (does not depend on time), this parameter is not used.

  \param[in] x[4]
  point in phase space, 4 coordinates: (l,L,g,G).

  \param[out] y[4]
  vectorfield at (t,x), 4 coordinates: d/dt(l,L,g,G).

  \param[in] params 
  pointer to the parameter of the system: the mass ratio "mu".

  \returns status code of the function (success/error).

  \retval GSL_SUCCESS success.
 */
int rtbp_del(double t, const double *x, double *y, void *params);

/**
  $l'$-component of Restricted Three Body Problem equations in Delaunay
  coordinates.

  \authors Marcel Guardia and Pau Roldan

  The function computes the $l'$-component of the vectorfield of the RTBP
  at a point named "x":
  \f[ \frac{d}{dt} l = L^{-3}+\mu\partial_L \Delta H_{circ}(x). \f]

  We normalize \f$\ell\f$ between 0 and 2\pi to compute the vectorfield.
  This is done so that function "eccentric" is more precise.

  \param[in] x
  point in phase space, 4 coordinates: (l,L,g,G).

  \param[out] dl 
  $l'$-component of the vectorfield at x.

  \param[in] params 
  pointer to the parameter of the system: the mass ratio "mu".

  \returns status code of the function (success/error).

  \retval GSL_SUCCESS success.

  \sa \ref rtbp_del
 */
int dot_l(const double *x, double *dl, void *params);

int dot_g(const double *x, double *dg, void *params);

/**
  Compute the function f0.

  The function computes the function f0, defined as
    \f[ f0(x) = 
    \frac{3-L^{-3}-\mu\partial_L \Delta H_{circ}(x)}
    {3(L^{-3}+\mu\partial_L \Delta H_{circ}(x))} = 
    \frac{3-l'(x)}{3 l'(x)}. \f]
  at a point named "x".

  \param[in] x		point in phase space, 4 coordinates: (l,L,g,G).
  \param[out] res	value of $f0$ (output)

  \param[in] params	
  pointer to the parameter of the system: the mass ratio "mu".

  \returns status code of the function (success/error)

  \retval GSL_SUCCESS	success.

  \sa \ref integrand_omega_in
 */
int f0(const double *x, double *res, void *params);

/**
  Compute the function f0 for the Stochastic paper.

  The function computes the function f0, defined as
    \f[ f0(x) = 
    \frac{\mu\partial_G \Delta H_{circ}(x)}
    {-1 + \mu\partial_G \Delta H_{circ}(x))} \f]
  at a point named "x".

  \param[in] x		point in phase space, 4 coordinates: (l,L,g,G).
  \param[out] res	value of $f0$ (output)

  \param[in] params	
  pointer to the parameter of the system: the mass ratio "mu".

  \returns status code of the function (success/error)

  \retval GSL_SUCCESS	success.

  \sa \ref integrand_omega_pm
 */
int f0_stoch(const double *x, double *res, void *params);

double re_DHell(const double *x, void *params);
double im_DHell(const double *x, void *params);
