// ====================================================================
// Hyperbolic Splitting (Eigenvalues/vects) Associated to a Fixed Point
// ====================================================================
// FILE:          $RCSfile: hyperdel.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2012-05-31 10:35:52 $
//
// FUNCTIONS
// =========
//
// hyper_del
// ---------
// Consider a 2D map $F$ and its derivative $DF$.
// Given a hyperbolic fixed point $p$ of $F$, compute its hyperbolic
// splitting, i.e. its unstable/stable eigenvalues and eigenvectors. 

#include <stdlib.h>     	// EXIT_SUCCESS, EXIT_FAILURE
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <prtbpdel.h>  		// section_t
#include <dprtbpdel_2d.h> 
#include <hyper.h> 		// ERR_DPRTBP_2D, ERR_NOTREAL

// name OF FUNCTION: hyper_del
//
// PURPOSE
// =======
// Consider a 2D map $F$ and its derivative $DF$.
// Given a hyperbolic fixed point $p$ of $F$, compute its hyperbolic
// splitting, i.e. its unstable/stable eigenvalues and eigenvectors. 
//
// PARAMETERS
// ==========
// mu
//    mass parameter for the RTBP
//
// sec
//    Poincare section
//
// H
//    energy value
//
// n
//    number of iterates of the Poincare map (n cuts with the
//    Poincare section).
// 
// p
//    Argument of the 2d derivative, 2 coordinates: p=(g,G). 
//
// eval
//    On exit, eval[0] holds the unstable eigenvalue, and eval[1] holds the
//    stable eigenvalue.
//
// evect
//    On exit, the first row evec[0,1] holds the unstable eigenvector, and
//    the second row evect[2,3] holds the stable eigenvector.
// 
// RETURN VALUE
// ============
// Returns a non-zero error code to indicate an error and 0 to indicate
// success.
// If an error is encountered, the function returns a non-zero value:
//
// ERR_DPRTBP_2D
//    Error computing derivative of 2D Poincare map of RTBP.
//
// ERR_NOTREAL
//    Eigenvalues are not real, so fixed point is not hyperbolic!
//
// NOTES
// =====
// Since the fixed point is hyperbolic, we assume the eigenvals/vects to be
// real.
//
// Since $DF$ is a 2x2 matrix, we could compute this simply by hand, see for
// instance
// http://www.math.harvard.edu/archive/21b_fall_04/exhibits/2dmatrices/index.html
//
// CALLS TO: dprtbp_del_2d

int hyper_del(double mu, section_t sec, double H, int n, double p[2], 
      double eval[2], double evec[4])
{
   double dp[4];        // derivative of 2D Poincare map

   gsl_vector_complex *evalc = gsl_vector_complex_alloc (2);
   gsl_matrix_complex *evecc = gsl_matrix_complex_alloc (2, 2);
   gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc (2);

   // auxiliary variables
   int status, i;
   gsl_complex eval_u, eval_s;			// unst/stable eigenvalue

   // Compute derivative of n-th iterate of 2D Poincare map, $DP^n(x)$.
   status=dprtbp_del_2d(mu,sec,H,n,p,dp);
   if(status)
   {
      fprintf(stderr, \
	    "hyper: error computing derivative of 2D Poincare map\n",n);
      return(ERR_DPRTBP_2D);
   }

   // Solve eigenvalue problem
   gsl_matrix_view m = gsl_matrix_view_array(dp,2,2);

   gsl_eigen_nonsymmv (&m.matrix, evalc, evecc, w);
   gsl_eigen_nonsymmv_free (w);

   // Only GSL_EIGEN_SORT_ABS_ASC and GSL_EIGEN_SORT_ABS_DESC are supported
   // due to the eigenvalues being complex.
   gsl_eigen_nonsymmv_sort (evalc, evecc, GSL_EIGEN_SORT_ABS_DESC);
   
   // Eigenvalues are assumed to be real.
   eval_u = gsl_vector_complex_get (evalc, 0);
   eval_s = gsl_vector_complex_get (evalc, 1);
   if(GSL_IMAG(eval_u)!=0 || GSL_IMAG(eval_s)!=0)
   {
      fprintf(stderr, "hyper: eigenvalues are not real!\n",n);
      return(ERR_NOTREAL);
   }
   eval[0] = GSL_REAL(eval_u);
   eval[1] = GSL_REAL(eval_s);

   // Eigenvectors are assumed to be real.
   gsl_vector_complex_view 
      evec_u = gsl_matrix_complex_column (evecc, 0);
   gsl_vector_complex_view 
      evec_s = gsl_matrix_complex_column (evecc, 1);
   evec[0] = GSL_REAL(gsl_vector_complex_get(&evec_u.vector, 0));
   evec[1] = GSL_REAL(gsl_vector_complex_get(&evec_u.vector, 1));
   evec[2] = GSL_REAL(gsl_vector_complex_get(&evec_s.vector, 0));
   evec[3] = GSL_REAL(gsl_vector_complex_get(&evec_s.vector, 1));

   gsl_vector_complex_free(evalc);
   gsl_matrix_complex_free(evecc);
     
   return 0;
}
