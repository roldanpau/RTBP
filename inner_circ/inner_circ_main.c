/*! \file
    \brief Inner Map of the Circular Problem: main prog.

    Compute the shift $T_0$ of the inner map of the circular problem.
    Moreover, this program also computes the partial shifts
    $\omega_{in}^{f,b}$ that we will need when computing the outer map.

    $Author: roldan $
    $Date: 2013-03-26 22:20:08 $
*/

#include <stdio.h>	// perror
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE

#include "inner_circ.h"	// inner_circ, omega_in_f, omega_in_b

/**
  Inner Map of the Circular Problem: main prog.

  This program computes the shift $T_0$ of the inner map of the circular
  problem, given by the integral
 
  \f[ 2\pi + T_0 = 
        \int_0^{6\pi} 
           \frac{1}{L^{-3}+\mu\partial_L \Delta H_{circ}(\lambda_I^1(s))} ds, \f]
 
  where \f$\gamma(s)\f$ is the periodic trajectory in the level of energy H,
  or equivalently
 
  \f[   T_0 = 
        \int_{0}^{6\pi} \frac{1}{\dot \ell (\lambda_I^1(s))} ds - 2\pi. \f]
 
  Morevover, this program also computes the partial shifts
  $\omega_{in}^{f,b}$ that we will need when computing the outer map:
 
  \f{eqnarray*}{
    \omega_{in}^f(I) &= \int_0^{4\pi} f0(\lambda_I^2(s)) ds. \\
    \omega_{in}^b(I) &= \int_0^{2\pi} f0(\lambda_I^1(s)) ds. 
  \f}
 
  (See notes "Kirkwood Gaps" by Marcel).
 
  It reads the following input from stdin:
 
  - mu	mass parameter for the RTBP
  - bInner (flag) 
  bInner == 0 means OUTER/LOWER splitting = "b",
  bInner == 1 means INNER/UPPER splitting = "f".
  
  Then it reads a sequence of input lines from stdin:
  - p[DIM] periodic point, 4 coordinates: (l,L,g,G).
  Notice that $p=p_2$ corresponds to INNER/UPPER separatrix, $p=p_1$ to
  OUTER/LOWER separatrix.
 
  For each input line, it writes the following output in stdout:
  - $w_{in}^*$ (where *=f,b depending on bInner flag).

  \remark
  The integrand is evaluated at an equispaced sequence of points 
  t\in(0,6\pi). 
  
 */

int main( )
{
   double mu;
   double p[DIM];	/* periodic point p=p2 or p1*/
   double T;		/* value of integral T_0*/
   double w_in;		/* value of integral \omega_{in}^* */

   //    bInner == 1 means INNER splitting = "f"
   //    bInner == 0 means OUTER splitting = "b", 
   int bInner;

   // aux vars
   double t;

   // Input mass parameter, bInner flag
   if(scanf("%le %d", &mu, &bInner)<2)
   {
      perror("main: error reading input");
      exit(EXIT_FAILURE);
   }

   // Input periodic point p=p2 or p1 from stdin.
   while(scanf("%le %le %le %le", p, p+1, p+2, p+3) == 4)
   {
      // Compute shift $T_0$ for inner map of circular problem.
      // T_0 has been computed through the period, so we now are only interested
      // in the partial shifts $\omega_in$.
      /*
      inner_circ(mu, p, &T);
      */

      if(bInner == 1)
      {
	 // Compute $\omega_in^f$.
	 omega_in_f(mu, p, &w_in);
      }
      else if (bInner == 0)
      {
	 // Compute $\omega_in^b$.
	 omega_in_b(mu, p, &w_in);
      }

      // Output result to stdout.
      //printf("%.15e\n", T);
      printf("%.15e\n", w_in);
   }
   exit(EXIT_SUCCESS);
}
