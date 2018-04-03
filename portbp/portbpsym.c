/*! \file portbpsym.c
    \brief Find a symmetric periodic orbit of the RTBP

    $Author: roldan $
    $Date: 2013-02-08 15:01:07 $
*/

#include <stdio.h>
#include <stdlib.h>     // EXIT_SUCCESS, EXIT_FAILURE
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_roots.h>
#include <prtbp_2d.h>
#include <dprtbp_2d.h>
#include <rtbp.h>	// DIM

/// desired accuracy for perpendicular condition PX=0
const double TOL_PERP=1.e-13;

const int MAXITER=100;		///< max number of solver iterations

const int ERR_SOLVER_STUCK=1;
const int ERR_MAXITER_PO=2;

/// parameters to function "dist"
struct dparams
{
   double mu;
   section_t sec;	///< Poincare section
   double H;
   int k;
};

double perp (double x, void *params);
double perp_df (double x, void *params);
void perp_fdf (double x, void *params, double *y, double *dy);

void
print_state_fdf (size_t iter, gsl_root_fdfsolver * s, double px)
{
 fprintf (stderr, "iter = %3lu x = % .10f, PX = % .4e\n", 
       iter, 
       gsl_root_fdfsolver_root(s),
       px);
}

/**
  Find a symmetric periodic orbit of the RTBP.

  This function refines a trajectory of the RTBP which is close to a
  symmetric periodic orbit, until a true periodic orbit is obtained.
  We look for an orbit that is symmetric wrt the horizontal axis, with
  initial condition corresponding to cutting the horizontal axis
  perpendicularly (x,px)=(x,0). After integrating for ONE FULL
  period, we request that the orbit cuts again the horizontal axis
  perpendiculary, i.e. with final condition (X,PX)=(X,0).

  Let "perp" be the function that tests the perpendicular condition: given
  init.cond. x, integrates for one full period and returns the final cond.
  PX.  We need to solve the function perp(x) := PX = 0.
  To refine the (approximate) root to this function, we use a 1D Newton
  method.

  \param[in] mu 	mass parameter of RTBP.
  \param[in] sec	Poincare section
  \param[in] H		energy value where we will look for periodic orbit
  \param[in] k		number of cuts with section

  \param[in,out] root
  On entry, approximate root to the "perp" function.
  On exit, it contains the true root (up to some tolerance TOL_PERP).

  \returns 
  Returns a non-zero error code to indicate an error and 0 to indicate
  success.
  If an error is encountered, "pt" is set to the best approximation to the
  fixed point found so far, and the function returns a non-zero value:
 
  ERR_SOLVER_STUCK
     Solver is unable to further improve the solution (up to desired
     accuracy TOL_FIXED_PT).

  ERR_MAXITER_PO
     Maximum number of solver iterations MAXITER reached without
     converging to a fixed point while looking for periodic orbit.
*/
// NOTES
// =====
// The (1-dimensional) root-finding method requires the solution of 1
// equation,
//    perp(x) = PX = 0,
// where (X,PX) is the iterate of (x,px): $(X,PX)=P(x,px)$.
// This is 1 equation in 1 variables, namely the variable $x$.
//
// We make use of both the perpendicular function "perp" and its derivative.  
//
// CALLS TO: perp, perp_df, perp_fdf

int portbpsym(double mu, section_t sec, double H, int k, double *root)
{
   const gsl_root_fdfsolver_type *T;
   gsl_root_fdfsolver *s;
  
   int status;
   int err = 0;		// return code of this function
   size_t iter = 0;
  
   struct dparams p = {mu,sec,H,k};	// params to perp function
   gsl_function_fdf f = {&perp, &perp_df, &perp_fdf, &p};
  
   double x,px;

   x=*root;

   T = gsl_root_fdfsolver_newton;	
   s = gsl_root_fdfsolver_alloc (T);
   gsl_root_fdfsolver_set (s, &f, x);
  
   do
     {
	iter++;
	status = gsl_root_fdfsolver_iterate (s);
	if (status)   /* check if solver is stuck */
	{
	   fprintf(stderr, "portbpsym: %s\n", gsl_strerror(status));
	   err=ERR_SOLVER_STUCK;
	   break;
	}

	// Set refined root
	x=gsl_root_fdfsolver_root(s);

	px=perp(x, &p);
	//print_state_fdf(iter,s,px);

	// Returns GSL_SUCCESS if the test |px|<TOL_PERP succeeds.
	// Returns GSL_CONTINUE otherwise.
	status = gsl_root_test_residual (px, TOL_PERP);
      }
    while (status == GSL_CONTINUE && iter < MAXITER);

   if(iter==MAXITER)
      err = ERR_MAXITER_PO;
   // else we have found the root
  
    // Set result
    *root=x;

    gsl_root_fdfsolver_free (s);

    return err;
}

// name OF FUNCTION: perp
// CREDIT:
// PURPOSE:
//    Computes the distance between a given point (x,px) and its iterate
//    under the Poincare map, (X,PX)=P(x,px).
//    Specifically, the function sets the distance vector "f" to
//       f=(f1,f2)=(X-x,PX-px).
//    When f(x,px)=(0,0), the root (x,px) is an initial condition for a
//    periodic trajectory.
//
// PARAMETERS:
// - ic approximate fixed point, 2 coordinates: (x,px).
// - params pointer to the parameters of the function: mu, sec, H
// - f value of the distance function, 2 coordinates: (f1,f2).
// 
// RETURN VALUE
// ============
// status code of the function (success/error). The function should return an
// appropriate error code to the GSL library if the function cannot be
// computed.
// 
// GSL_SUCCESS 
//    success.
//
// ERR_POINCARE_MAP
//    Error computing the 2D Poincare map. The function cannot be computed.
//
// NOTES
// =====
// The 2D poincare map $P^{k}(x,px)$ is computed by calling the external
// function "prtbp_2d".
//
// CALLS TO: prtbp_2d
//
// CALLED FROM: portbpsym

double perp (double x, void *params)
{
   double mu=((struct dparams *)params)->mu;
   section_t sec=((struct dparams *)params)->sec;
   double H=((struct dparams *)params)->H;
   int k=((struct dparams *)params)->k;
   double pt[2];

   // this is not used, but it's needed in call to prtbp_2d
   double ti;	// integration time

   pt[0] = x;
   pt[1] = 0;	// px
   
   // Compute the 2D poincare map $P^(pt)$
   if(prtbp_2d(mu,sec,H,k,pt,&ti))
   {
      fprintf(stderr, "perp: error computing 2D poincare map\n");
      exit(EXIT_FAILURE);
   }
   
   // Return PX
   return(pt[1]);
}

// name OF FUNCTION: dist_df
// CREDIT:
// PURPOSE:
//    Computes Jacobian of the distance function "dist(x,px)" defined above.
//    Specifically, the function sets the (2x2) Jacobian matrix "J" to
//       J = DP - I,
//    where DP is the derivative of the Poincare map and I is the identity.
//
// PARAMETERS:
// - ic approximate fixed point, 2 coordinates: (x,px).
// - params pointer to the parameters of the function: mu, sec, H
// - J value of the 2x2 Jacobian matrix.
// 
// RETURN VALUE
// ============
// status code of the function (success/error). The function should return an
// appropriate error code to the GSL library if the Jacobian cannot be
// computed.
// 
// GSL_SUCCESS 
//    success.
//
// ERR_DPOINCARE_MAP
//    Error computing the derivative of the 2D Poincare map. The function
//    cannot be computed.
//
// NOTES
// =====
// The derivative DP of the 2D poincare map $P(x,px)$ is computed by calling
// the external function "dprtbp_2d".
//
// CALLS TO: dprtbp_2d
//
// CALLED FROM: portbp

double perp_df (double x, void *params)
{
   double mu=((struct dparams *)params)->mu;
   section_t sec=((struct dparams *)params)->sec;
   double H=((struct dparams *)params)->H;
   int k=((struct dparams *)params)->k;
   double pt[2];
   double dp[4];	// derivative of Poincare map

   pt[0] = x;
   pt[1] = 0; 	// px
   
   // Compute the derivative of 2D poincare map $P(pt)$
   if(dprtbp_2d(mu,sec,H,k,pt,dp))
   {
      fprintf(stderr, "perp: error computing derivative of 2D poincare map\n");
      exit(EXIT_FAILURE);
   }
   
   // Return \partial PX / \partial x
   return(dp[2]);
}

void perp_fdf (double x, void *params, double *y, double *dy)
{
   *y = perp(x,params);
   *dy = perp_df(x,params);
}
