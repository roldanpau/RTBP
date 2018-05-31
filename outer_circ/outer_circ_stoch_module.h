/*! \file outer_circ_stoch_module.h
    \brief Outer Map of the Circular Problem for Stochastic paper
*/

#ifndef OUTER_CIRC_STOCH_MODULE_H_INCLUDED
#define OUTER_CIRC_STOCH_MODULE_H_INCLUDED

#include <rtbp.h>   // DIM
#include <section.h>

/**
  Given an energy level \f$H\f$, compute \f$\omega_-^j(H)\f$.

  Given an energy level \f$H\f$, compute \f$\omega_-^j(H)\f$, defined as
 
  \f[ \omega_-^j(I) = 
        lim_{N\to -\infty} 
           \left( \int_0^{2\pi N} f0(\gamma^*(s)) ds + NT_0(I) \right),
  \f]
 
  where \f$\gamma^*(s)\f$ is the homoclinic trajectory of the reduced flow that
  starts at the primary homoclinic point z.
 
  Equivalently, we compute it as
  
  \f[\omega_-^j(I) = 
        lim_{N\to \infty} 
           \left( \int_{2\pi N}^0 f0(\gamma^*(s)) ds - NT_0(I) \right),
  \f]
 
  where \f$\gamma^*(s)\f$ is the homoclinic trajectory that starts at the point 
  \f$ z^u = \mathcal{P}^{-N}(z) \f$. Notice that z^u is in the UNSTABLE
  manifold of the reduced flow.
 
  Equivalently, we split the integral from 2\pi N to 0 into N parts of size
  2\pi:
 
  \f[\omega_-^j(I) = 
        lim_{N\to \infty} 
           \left( \sum_{i=N,1} 
              (\int_{2\pi}^{0} f0(\gamma_i(s)) ds - T_0(I)) \right),
  \f]
 
  where \f$\gamma_i(s)\f$ is the homoclinic trajectory that starts at the point
  P^{N-i}(z^u) = P^{-i}(z). 
 
  This is computed using numerical integration.
 
  \param[in] mu 	mass parameter for the RTBP
  \param[in] sec    Poincare section: sec={SECg,SECg2}

  \param[in] x    [DIM]	x=(l,L,g,G),     point z^u, on the section g=0.
  \param[in] x_car[DIM]	x=(x,y,p_x,p_y), point z^u, on the section g=0.

  \param[in] N
     number of iterates of the poincare map (length of integration).
     N must be a POSITIVE integer.

  \param[in] T0	shift of inner map

  \param[out] omega
     On return of this function, omega contains the value of the function
     \f$\omega_-^j\f$.
  
  \returns
  a non-zero error code to indicate an error and 0 to indicate success.
 */

int omega_neg_stoch(double mu, section_t sec, double x[DIM], double x_car[DIM], 
        int N, double T0, double *omega);

/**
  Given an energy level \f$H\f$, compute \f$\omega_+^j(H)\f$.

  Given an energy level \f$H\f$, compute \f$\omega_+^j(H)\f$, defined as
 
  \f[ \omega_+^j(I) = 
        lim_{N\to \infty} 
           \int_0^{2\pi N} f0(\gamma^*(s)) ds + NT_0(I),
  \f]
 
  where \f$\gamma^*(s)\f$ is the homoclinic trajectory of the reduced flow that
  starts at the primary homoclinic point z.
 
  Equivalently, we compute it as
  
  \f[\omega_+^j(I) = 
        lim_{N\to +\infty} 
           \int_{-2\pi N}^0 f0(\gamma^*(s)) ds + NT_0(I),
  \f]
 
  where \f$\gamma^*(s)\f$ is the homoclinic trajectory that starts at the point 
  \f$ z^s = \mathcal{P}^{N}(z) \f$. Notice that z^s is in the STABLE manifold of
  the reduced flow.
 
  Equivalently, we split the integral from -2\pi N to 0 into N parts of size
  2\pi:
 
  \f[\omega_+^j(I) = 
        lim_{N\to +\infty} 
           \sum_{i=N,1} 
              (\int_{-2\pi}^{0} f0(\gamma_i(s)) ds + T_0(I)),
  \f]
 
  where \f$\gamma_i(s)\f$ is the homoclinic trajectory that starts at the point
  P^{-(N-i)}(z^s) = P^{i}(z). 
 
  This is computed using numerical integration.
 
  \param[in] mu 	mass parameter for the RTBP
  \param[in] sec    Poincare section: sec={SECg,SECg2}

  \param[in] x    [DIM]	x=(l,L,g,G),     point z^s, on the section g=0.
  \param[in] x_car[DIM]	x=(x,y,p_x,p_y), point z^s, on the section g=0.

  \param[in] N
     number of iterates of the poincare map (length of integration).
     N must be positive since we are computing \f$\omega_+^j\f$.

  \param[in] T0	shift of inner map

  \param[out] omega
     On return of this function, omega contains the value of the function
     \f$\omega_+^j\f$.
  
  \returns
  a non-zero error code to indicate an error and 0 to indicate success.
 */

int omega_pos_stoch(double mu, section_t sec, double x[DIM], double x_car[DIM],
		int N, double T0, double *omega);

#endif // OUTER_CIRC_STOCH_MODULE_H_INCLUDED

