/*! \file portbp.c
    \brief Find a periodic orbit of the RTBP

    $Author: roldan $
    $Date: 2013-02-08 15:01:07 $
*/

#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

#include <prtbp.h>		// section_t
#include <prtbp_2d.h>
#include <dprtbp_2d.h>
#include <rtbp.h>	// DIM

const double TOL_FIXED_PT=1.e-13;	///< desired accuracy for fixed point
const int MAXITER=100;			///< max number of solver iterations

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

int dist (const gsl_vector * ic, void *params, gsl_vector * f);
int dist_df (const gsl_vector * ic, void *params, gsl_matrix * J);
int dist_fdf (const gsl_vector * ic, void *params, gsl_vector * f, 
      gsl_matrix * J);

void
print_state (size_t iter, gsl_multiroot_fdfsolver * s)
{
 printf ("iter = %3u p = % .10f % .10f "
	 "dist(p) = % .3e % .3e\n",
	 iter,
	 gsl_vector_get (s->x, 0),
	 gsl_vector_get (s->x, 1),
	 gsl_vector_get (s->f, 0),
	 gsl_vector_get (s->f, 1));
}

/**
  Find a periodic orbit of the RTBP.

  This function refines a trajectory of the RTBP which is close to a periodic
  orbit, until a true periodic orbit is obtained.
  We look for a periodic orbit in the given energy manifold $H=H_0$.
  Consider the Poincare section "sec" in the autonomous RTBP.
  We look for the periodic periodic orbit as a fixed point of the $k$-th
  iterate $P^{k}(x)$ of the Poincare first return map to this section.
  To refine the (approximate) fixed point, we use a Newton-like method.

  \param[in] mu 	mass parameter of RTBP.
  \param[in] sec	Poincare section
  \param[in] H		energy value where we will look for periodic orbit
  \param[in] k		number of cuts with section

  \param[in,out] pt
  On entry, approximate fixed point of the $k$-th iterate $P^{k}(x)$. 
  On exit, it contains the true fixed point.
    
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
// The (2-dimensional) root-finding method requires the simultaneous solution
// of 2 equations,
//    dist_1(x,px) = X - x = 0
//    dist_2(x,px) = PX - px = 0,
// where (X,PX) is the iterate of (x,px): $(X,PX)=P(x,px)$.
// These are 2 equations in 2 variables: x,px.
//
// We make use of both the distance function "dist" and its derivative.  
//
// A standard multidimensional solver from the GSL Scientific Library is
// used, namely gsl_multiroot_fdfsolver_hybridj. This is a modified version
// of Powell's Hybrid method (see the GSL manual for details) without scaling.
// It is important to note that we have tried other numerical solvers, such
// as gsl_multiroot_fsolver_hybrids and gsl_multiroot_fsolver_dnewton, and
// they did not work as well (failed to converge for some energy values).

int portbp(double mu, section_t sec, double H, int k, double pt[2])
{
   const gsl_multiroot_fdfsolver_type *T;
   gsl_multiroot_fdfsolver *s;
  
   int status;
   int err = 0;		// return code of this function
   size_t iter = 0;
  
   const size_t n = 2;
   struct dparams p = {mu,sec,H,k};	// params to dist function
   gsl_multiroot_function_fdf f = {&dist, &dist_df, &dist_fdf, n, &p};
  
   gsl_vector *x = gsl_vector_alloc (n);

   gsl_vector_set (x, 0, pt[0]);
   gsl_vector_set (x, 1, pt[1]);
  
   // hybrid method without internal scaling is the one that works best
   T = gsl_multiroot_fdfsolver_hybridj;	
   s = gsl_multiroot_fdfsolver_alloc (T, n);
   gsl_multiroot_fdfsolver_set (s, &f, x);
  
   print_state (iter, s);
  
   do
     {
	iter++;
	status = gsl_multiroot_fdfsolver_iterate (s);
  
	print_state (iter, s);
  
	if (status)   /* check if solver is stuck */
	{
	   fprintf(stderr, "portbp: %s\n", gsl_strerror(status));
	   err=ERR_SOLVER_STUCK;
	   break;
	}
	status = gsl_multiroot_test_residual (s->f, TOL_FIXED_PT);
      }
    while (status == GSL_CONTINUE && iter < MAXITER);

   if(iter==MAXITER)
      err = ERR_MAXITER_PO;
   // else we have found the fixed point
  
    // Set refined fixed point
    pt[0]= gsl_vector_get (s->x, 0);
    pt[1]= gsl_vector_get (s->x, 1);
  
    gsl_multiroot_fdfsolver_free (s);
    gsl_vector_free (x);

    return err;
}

// name OF FUNCTION: dist
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
// CALLED FROM: portbp

const int ERR_POINCARE_MAP=1;

int dist (const gsl_vector * ic, void *params, gsl_vector * f)
{
   double mu=((struct dparams *)params)->mu;
   section_t sec=((struct dparams *)params)->sec;
   double H=((struct dparams *)params)->H;
   int k=((struct dparams *)params)->k;
   double x_init[2];	// initial condition, x_init=(x,px)
   double pt[2];

   // this is not used, but it's needed in call to prtbp_2d
   double ti;	// integration time

   x_init[0] = pt[0] = gsl_vector_get(ic,0);
   x_init[1] = pt[1] = gsl_vector_get(ic,1);
   
   printf("pt=%e %e\n", pt[0], pt[1]);
   // Compute the 2D poincare map $P^(pt)$
   if(prtbp_2d(mu,sec,H,k,pt,&ti))
   {
      fprintf(stderr, "dist: error computing 2D poincare map\n");
      return(ERR_POINCARE_MAP);
   }
   
   // Set the distance vector "f" to f=(f1,f2)=(X-x,PX-px).
   gsl_vector_set(f,0,pt[0]-x_init[0]);
   gsl_vector_set(f,1,pt[1]-x_init[1]);
   return(GSL_SUCCESS);
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

const int ERR_DPOINCARE_MAP=1;

int dist_df (const gsl_vector * ic, void *params, gsl_matrix * J)
{
   double mu=((struct dparams *)params)->mu;
   section_t sec=((struct dparams *)params)->sec;
   double H=((struct dparams *)params)->H;
   int k=((struct dparams *)params)->k;
   double x_init[2];	// initial condition, x_init=(x,px)
   double pt[2];
   double dp[4];	// derivative of Poincare map

   x_init[0] = pt[0] = gsl_vector_get(ic,0);
   x_init[1] = pt[1] = gsl_vector_get(ic,1);
   
   // Compute the derivative of 2D poincare map $P(pt)$
   if(dprtbp_2d(mu,sec,H,k,pt,dp))
   {
      fprintf(stderr, "dist: error computing derivative of 2D poincare map\n");
      return(ERR_DPOINCARE_MAP);
   }
   
   // Set the Jacobian matrix "J" to J=dp-I
   gsl_matrix_set(J,0,0,dp[0]-1.0); gsl_matrix_set(J,0,1,dp[1]);
   gsl_matrix_set(J,1,0,dp[2]);     gsl_matrix_set(J,1,1,dp[3]-1.0);
   return(GSL_SUCCESS);
}

int dist_fdf (const gsl_vector * ic, void *params, gsl_vector * f, 
      gsl_matrix * J)
{
   dist(ic,params,f);
   dist_df(ic,params,J);
   return GSL_SUCCESS;
}
