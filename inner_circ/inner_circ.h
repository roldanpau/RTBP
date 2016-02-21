/*! \file
    \brief Inner Map of the Circular Problem.

    $Author: roldan $
    $Date: 2013-03-26 22:20:08 $
*/

#include <rtbp.h>	// DIM

/// Parameters to the \ref integrand_omega_in function.
struct iparams_omega_in
{
   double mu;
   double x[DIM];
};

/**
  Integrand of \f$ \omega_{in} \f$ integral.

  Given a trajectory \f$ \lambda(s) \f$ of the reduced circular RTBP, and an
  integration time 's' in the trajectory, this function computes the
  integrand
    \f[ f0(\lambda(s)) = 
    \frac{3-L^{-3}-\mu(\partial_L \Delta H_{circ}) \circ \lambda(s)}
    {3(L^{-3}+\mu(\partial_L \Delta H_{circ}) \circ \lambda(s))} \f]
  evaluated at \f$\lambda(s)\f$.

  \param[in] s	integration time in the periodic trajectory.
  \param[in] mu	mass parameter for the RTBP
  \param[in] x[DIM]	periodic point, 4 coordinates: (l,L,g,G).

  \returns	the integrand evaluated at the point \f$\lambda(s)\f$.

  \remark
  The parameter s may be positive or negative.

  \remark
  We will also use this as the integrand in \ref omega_pos and \ref omega_neg.

  \sa \ref f0
 */

double integrand_omega_in(double s, void *params);

/**
  Compute \f$\omega_{in}^f(I)\f$

  Given an energy level $H=-I$, compute \f$\omega_{in}^f(I)\f$, defined as
  \f[   \omega_{in}^f(I) = 
        \int_0^{4\pi} f0(\lambda_I^2(s)) ds. 
  \f]
  This is computed using numerical integration.

  \param[in] mu	mass parameter for the RTBP
  \param[in] x[DIM]
  x=(l=0,L,g,G), periodic point of period 3. This point corresponds to the
  unique periodic orbit for the circular problem with initial condition at
  \f$ \tilde\Lambda_0^2 \f$ and energy $H=-I$. 
  Note that x is on the section l=0.

  \param[out] omega
  On return of this function, omega contains the value of the integral.

  \returns
  a non-zero error code to indicate an error and 0 to indicate success.
  
 */
int omega_in_f(double mu, double x[DIM], double *omega);
int omega_in_b(double mu, double x[DIM], double *omega);

/**
  Inner Map of the Circular Problem.

  For the circular problem, the invariant cylinder is foliated by periodic
  orbits.
  It is assumed that we work with a 3:1 resonant periodic orbit that is
  symmetric.
  According to Marcel's notes, the inner map for the circular problem is given
  by
  
  \f[ F_0^{in}\colon\ (I,t)\to (I, t+ T_0(I)). \f]
 
  Here, $T_0(I)$ is given by the integral
  
  \f[ 2\pi + T_0 = 
        \int_0^{6\pi} 
           \frac{1}{L^{-3}+\mu\partial_L \Delta H_{circ}(\gamma(s))} ds, \f]
 
  where \f$\gamma(s)\f$ is the periodic trajectory in the level of energy H.
 
  Given an energy level $H$, this function computes $T_0(H)$.
 
  \param[in] mu
     mass parameter for the RTBP

  \param[in] x
     x=(l=0,L,g,G), periodic point of period 3, on the section l=0.

  \param[out] T
     On return of this function, T contains the shift T_0(H).
  
  \returns
  a non-zero error code to indicate an error and 0 to indicate success.

  \remark
  We recall that \f$2\pi+T_0\f$ in the original time is just the period of the
  periodic orbit, which we have already computed using program "porbits". 
  As a test, here we compute T_0 using the integral above.
  
 */
int inner_circ(double mu, double x[DIM], double *T);
