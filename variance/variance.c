/*! \file
    \brief Compute the first order of the variance $\sigma_0^2$: main prog.

	See the formula for the first order of the variance in Ansatz 2 of the
	second paper Stochastic3BP/Paper2NHIL3BP.
*/

#include <stdio.h>	// perror
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE
#include <string.h> // strcmp
#include <assert.h>
#include <math.h>   // cos, sin
#include <complex.h>	// Standard Library of Complex Numbers

#include "variance_module.h"	// Mean

const int N_THETA = 30;	// number of samples for variable \theta\in[0,2\pi).

/**
  Compute the first order of the variance $\sigma_0^2$: main prog.

  It reads the following input from stdin:

  A sequence of lines:

  - alpha_1, alpha_2
     phase shift of homoclinic channels 1 and 2
  - Re B_1, Im B_1, Re B_2, Im B_2
     Melnikov functions of homoclinic channels 1 and 2 (real and imaginary
	 parts).

  For each input line, it outputs result to stdout:
  - variance
     first order of the variance.
 
 */
 
int main( )
{
	double H;

	double alpha_neg_1, alpha_pos_2;
	double alpha1, alpha2;

	double ReB1, ImB1;
	double ReB2, ImB2;

	double complex B1;
	double complex B2;

	double beta1, beta2;	// \beta_i^0


	// auxiliary variables
	double mod1;
	double mod2;
	double theta;
	double complex exp1, exp2;	// exp{i \beta_i^0}

	// Input alpha_neg_1, alpha_pos_2, Re(B_1), Im(B_1), Re(B_2), Im(B_2), from
	// stdin.
	while(scanf("%le %le %le %le %le %le %le", &H, &alpha_neg_1, &alpha_pos_2,
				&ReB1, &ImB1, &ReB2, &ImB2) == 7)
	{
		alpha1 = -2*alpha_neg_1;
		alpha2 = +2*alpha_pos_2;

		for(int i=0; i<=N_THETA; i++)
		{
			theta = i*(2*M_PI/N_THETA);

			beta1 = theta + alpha1;
			beta2 = theta + alpha2;

			exp1 = cos(beta1) + I*sin(beta1);
			exp2 = cos(beta2) + I*sin(beta2);

			B1 = ReB1 + I*ImB1;
			B2 = ReB2 + I*ImB2;

			mod1 = cabs(B1 - Mean(B1,B2) * (1 - exp1) / (1 - Mean(exp1, exp2)));
			mod2 = cabs(B2 - Mean(B1,B2) * (1 - exp2) / (1 - Mean(exp1, exp2)));

			printf("%f %f %.2e\n", H, theta, mod1*mod1 + mod2*mod2); 
		}
		printf("\n");
	}

	exit(EXIT_SUCCESS);
}
