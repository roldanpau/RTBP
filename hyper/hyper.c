/*! \file
    \brief Hyperbolic Splitting (Eigenvalues/vects) Associated to a Fixed
    Point

    $Author: roldan $
    $Date: 2013-03-11 11:22:18 $
*/

#include <stdio.h>  	// fprintf
#include <math.h>  	// pow,sqrt
#include <assert.h> 
#include <section.h>  	// section_t
#include <dprtbp_2d.h>  // dprtbp2d

const int ERR_DPRTBP_2D=1;
const int ERR_NOTREAL=2;

int hyper(double mu, section_t sec, double H, int n, double p[2], 
      double eval[2], double evec[4])
{
   double dp[4];        // derivative of 2D Poincare map
   double L1,L2;	// eigenvalues

   // auxiliary variables
   double T;		// trace of DP matrix
   double D;		// determinant of DP matrix
   double disc;		// discriminant
   double a,b,c,d;	// entries in DP matrix
   int status;
   double norm;		// norm of eigenvector

   // Compute derivative of n-th iterate of 2D Poincare map, $DP^n(x)$.
   status=dprtbp_2d(mu,sec,H,n,p,dp);
   if(status)
   {
      fprintf(stderr, \
	    "hyper: error computing derivative of 2D Poincare map\n");
      return(ERR_DPRTBP_2D);
   }

   a=dp[0];	b=dp[1];
   c=dp[2];	d=dp[3];

   T=a+d;	// trace
   D=a*d-b*c;	// determinant

   // Eigenvalues
   disc = pow(T,2)/4-D;
   //fprintf(stderr, "a=%.8e b=%.8e c=%.8e d=%.8e disc=%.8e\n", a,b,c,d,disc);
   if(disc<0)
   {
      fprintf(stderr, "hyper: eigenvalues are not real!\n");
      return(ERR_NOTREAL);
   }
   L1 = T/2 + sqrt(disc);
   L2 = T/2 - sqrt(disc);

   // Notice that L1,L2>0 and L1>L2. In fact, L1>1, L2<1.
   assert(L1>1 && 0<L2 && L2<1);
   eval[0] = L1;	// unstable eigenvalue
   eval[1] = L2;	// stable eigenvalue

   // Eigenvectors
   if(c!=0)
   {
      // unst eigenvector
      evec[0] = L1-d;
      evec[1] = c;

      // st eigenvector
      evec[2] = L2-d;
      evec[3] = c;
   }
   else if(b!=0)	// c=0 and b!=0
   {
      // unst eigenvector
      evec[0] = b;
      evec[1] = L1-a;

      // st eigenvector
      evec[2] = b;
      evec[3] = L2-a;
   }
   else			// both b and c are 0
   {
      // unst eigenvector
      evec[0] = 1;
      evec[1] = 0;

      // st eigenvector
      evec[2] = 0;
      evec[3] = 1;
   }

   // Normalise eigenvectors to unit magnitude.
   norm = sqrt(pow(evec[0],2)+pow(evec[1],2));
   evec[0] /= norm;
   evec[1] /= norm;

   norm = sqrt(pow(evec[2],2)+pow(evec[3],2));
   evec[2] /= norm;
   evec[3] /= norm;

   // we assume that eigenvector $v=(x,p_x)$ points "to the right", i.e. we
   // assume that the first component of $v$ is $x>0$ (see approxint_unst).
   if(evec[0]<0)
   {
      evec[0] = -evec[0];
      evec[1] = -evec[1];
   }
   if(evec[2]<0)
   {
      evec[2] = -evec[2];
      evec[3] = -evec[3];
   }
   return 0;
}
