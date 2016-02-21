/*! \file Ldeviation.c
    \brief Deviation of semi-major axis L wrt resonant one
    $Author: roldan $
    $Date: 2013-02-08 14:46:48 $
*/

#include <stdio.h>	// fprintf
#include <math.h>	// fmod, fabs, cbrt
#include <gsl/gsl_errno.h>      // GSL_SUCCESS
#include <frtbpdel.h>
#include <rtbp.h>	// DIM

// M_PI has dropped from the ISO C std, so we define it ourselves
#ifndef M_PI
#define M_PI 3.14159265358979323844
#endif 

/**
Given a periodic orbit \f$\gamma_J(t)\f$ of the RTBP in rotating Delaunay
coordinates $(l,L,g,G)$, compute the maximum deviation of the $L$ coordinate
(sqrt of semi-major axis) from that of the true resonant periodic orbit, i.e.
compute
   \f[ max{|L(t)-L^*|} \f]
for all \f$t \in [0,T_J]\f$, where $L^* = 3^{-1/3}$ for the 3:1 resonance.

\param[in] mu mass parameter for the RTBP.
\param[in,out] x
   Initial point, 4 coordinates: (l, L, g, G). 
   On return of the this function, this point is modified!
\param[out] Ldev
   On return, it holds the maximum deviation of $L$ from $L^*$.

\return Returns a non-zero error code to indicate an error and 0 to indicate
success.
If an integration error is encountered, the function returns a non-zero
value.
*/
int Ldeviation(double mu, double x[DIM], double *Ldev)
{
   const double Lstar = pow(3,-1.0/3.0);

   double x_pre[DIM];   /* previous value of point x */
   int status;
   int i,n;


   // auxiliary variables
   double n1,n2;
   double L;		// angular momentum
   double Ldev_new;

   double t=0;
   (*Ldev) = 0.0;
   while(t<2*M_PI)
   {
	 //printf("x = %le %le %le %le\n",x[0], x[1], x[2], x[3]);

         // Save previous value of point "x"
         for(i=0;i<DIM;i++)
            x_pre[i]=x[i];

	 // Integrate for a "short" time t1=0.1, short enough so that we can
	 // detect crossing of Poincare section.

	 // WARNING! Before we used t1=1 as a "short" time, but sometime this
	 // was too long...
	 status = frtbp_del(mu,0.001,x);
	 t+=0.001;
	 if (status != GSL_SUCCESS)
	 {
	    fprintf(stderr, "Ldeviation: error integrating trajectory\n");
	    return(1);
	 }

	 // Update Ldev
	 L=x[1];
	 //Ldev_new = fabs(L-Lstar);
	 Ldev_new = L-Lstar;
	 if(Ldev_new > (*Ldev))
	    (*Ldev) = Ldev_new;
	 
   }
   return(0);
}
