/*! \file intersec.c
    \brief Intersection of Invariant Manifolds

    $Author: roldan $
    $Date: 2013-03-11 11:27:15 $
*/

#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_roots.h>

#include <rtbp.h>	// DIM
#include <hinv.h>
#include <prtbp_nl.h>	// prtbp_nl, prtbp_nl_inv
#include <prtbp_nl_2d_module.h>	// prtbp_nl_2d, prtbp_nl_2d_inv

/// Tolerance (precision) for bisection method 
const double BISECT_TOL=1.e-15;

// Parameters to distance_f_unst and distance_f_st functions.
struct dparams
{
   double mu;		// mass parameter
   double H;		// energy value
   double p[2];		// fixed point
   double v[2];		// unstable/stable vector
   double n;		// num. of iteration in the unstable/stable dir.
   double l;		// axis line
};

double
distance_f_unst (double h, void *params);
double
distance_f_st (double h, void *params);

int
print_state (size_t iter, gsl_root_fsolver * s)
{
  double r = gsl_root_fsolver_root (s);
  double x_lo = gsl_root_fsolver_x_lower (s);
  double x_hi = gsl_root_fsolver_x_upper (s);

 fprintf (stderr, "iter = %3lu, [%.7f, %.7f], root = %.7f, err(est) = %.15f\n",
	 iter,
	 x_lo, 
	 x_hi,
	 r,
	 x_hi-x_lo);
}

/**
  Intersection of unstable invariant manifold with symmetry line.

  Consider the 2D map \f$\mathcal{P}: \Sigma_- \to \Sigma_-\f$, which is
  assumed to be reversible with respect to the symmetry line $p_x=0$.
  Let $p$ be a hyperbolic fixed point for \f$\mathcal{P}\f$.
  For definiteness, we assume that $p$ is located above the $p_x=0$ axis.

  Assume that \f$\lambda\f$ is the unstable eigenvalue, with \f$\lambda>1\f$.
  Let $v$ be the unstable eigenvector for the eigenvalue \f$\lambda\f$. 
  For definiteness, we assume that $v=(x,p_x)$ points "to the right", i.e. we
  assume that the first component of $v$ is $x>0$, but we could use the other
  branch of the manifold.
  Let $W^u(p)$ be the unstable manifold of $p$.
  Let $p_x=l$ be a line parallel to the $x$ axis.

  This function computes the "first" intersection point $z=(x,l)$ of the
  unstable manifold with the line $p_x=l$ as we grow the manifold from the
  fixed point.
  It does so by refining the approximate homoclinic point obtained in \ref
  approxint_unst.

  We consider the $n$-th iteration under the iterated Poincare map of the
  fundamental segment between the two points $p+h v$ and $P(p+h v)$.
  Let $p_u = p+h_u v$ be a point in the unstable segment.
  We are provided (by the program approxint) with an interval $(h_1,h_2)$
  that is guaranteed to contain the root $h^*$ such that 
     \f[ p_x(P^n(p+h^* v))=l. \f]
  Then we use a standard numerical method (bisection-like 1-dimensional root
  finding) to find $h^*$.
  The function that we solve (find a zero of) is
     distance(h_u) = l - p_x(P^n(p+h_u v_u)).

  \remark
  We ask for precision of BISECT_TOL in the homoclinc point.
  
  \param[in] mu         mass parameter for the RTBP
  \param[in] H          energy value
  \param[in] p          fixed point $p=(x,p_x)$
  \param[in] v          eigenvector associated to unstable direction
  \param[in] lambda     eigenvalue associated to unstable direction

  \param[in] n
  number of desired iterations by the Poincare map in the unstable direction.

  \param[in] h1,h2
  Endpoints of the unstable segment u_i bracketing the approximate root
  $p_u=p + h_u v_u$ with $h_u \in (h_1,h_2)$ (initial interval for bisection
  method). This is furnished by function \ref approxint_unst.

  \param[in] l          line $p_x=l$ parallel to the $x$ axis.

  \param[out] h
  On exit, it contains the true root found by numerical method.

  \param[out] p_u[2]
  On exit, it contains point in the unstable segment such that
  p_x(P(p_u)) = l.

  \param[out] tu
  On exit, it contains integration time to reach z from p_u.
  
  \param[out] z[DIM]	
  On exit, it contains homoclinic point z = P^n(p_u).

  \returns 
  a non-zero error code to indicate an error and 0 to indicate
  success.

  \retval 1 	Problems computing the Poincare iterates.
  \retval 2 	Bisection procedure did not converge
*/

// NOTES
// =====
// For the moment, we work with the 3:1 resonant family of periodic orbits.

int intersec_unst(double mu, double H, double p[2], double v[2], 
      double lambda, int n, double h1, double h2, double l,
      double *h, double p_u[2], double *t, double z[DIM])
{
   const gsl_root_fsolver_type *T;
   gsl_root_fsolver *s;

   struct dparams params;
   params.mu = mu;
   params.H = H;
   params.p[0] = p[0];
   params.p[1] = p[1];
   params.v[0] = v[0];
   params.v[1] = v[1];
   params.n = n;
   params.l = l;
   gsl_function f = {&distance_f_unst, &params};

   // auxiliary variables
   double htmp;

   // For some reason, bisection method complains that interval [h1,h2] does
   // not straddle 0, so we try to enlarge it a little bit.
   //h1 *= 0.999;
   //h2 /= 0.999;

   if(h2<h1)
   {
	   htmp = h1;
	   h1 = h2;
	   h2 = htmp;
   }

   T = gsl_root_fsolver_brent;
   s = gsl_root_fsolver_alloc (T);
   gsl_root_fsolver_set (s, &f, h1, h2);

   // auxiliary vars
   int status;
   size_t iter = 0;
   double x_lo, x_hi;

   // Find a root of the distance function, i.e. an intersection point of
   // the manifolds (using a bisection method).
   //print_state (iter, s);

    do
      {
	iter++;
	status = gsl_root_fsolver_iterate (s);
	if (status)   /* check if solver is stuck */
	  break;
  
	x_lo = gsl_root_fsolver_x_lower(s);
	x_hi = gsl_root_fsolver_x_upper(s);

	// epsabs=BISECT_TOL, epsrel=0
	status =
	  gsl_root_test_interval (x_lo, x_hi, BISECT_TOL, 0);
        //print_state (iter, s);
      }
    while (status == GSL_CONTINUE && iter < 1000);

    //fprintf (stderr, "status = %s\n", gsl_strerror (status));

    // the root is:
    *h = gsl_root_fsolver_root(s);

    // If bisection did not converge, warn calling function.
    // In this case, the root is updated to the closest zero, 
    // but we don't compute further results, which would be unaccurate (ps,
    // pu, z)
    if(status != GSL_SUCCESS)
       return(2);

   // Compute the following:
   // - point p_u

    p_u[0] = p[0] + (*h) * v[0];
    p_u[1] = p[1] + (*h) * v[1];

   // - intersection point z = P(p_u),
   // - t: integration time to reach homoclinic point $z$.

   z[0] = p_u[0];	// x
   z[1] = 0;		// y
   z[2] = p_u[1];	// px
   status=hinv(mu,SEC2,H,z);
   if(status)
   {
      fprintf(stderr, "intersec: error lifting point\n");
      return(1);
   }

   status=prtbp_nl(mu,SEC2,4*n,z,t); 	// $z = P^n(z)$
   if(status)
      {
	 fprintf(stderr, "intersec: error computing intersection point\n");
	 return(1);
      }
   return(0);
}

/**
  Intersection of stable invariant manifold with symmetry line.

  Exactly as \ref intersec_unst.
  */

int intersec_st(double mu, double H, double p[2], double v[2], 
      double lambda, int n, double h1, double h2, double l,
      double *h, double p_s[2], double *t, double z[DIM])
{
   const gsl_root_fsolver_type *T;
   gsl_root_fsolver *s;

   struct dparams params;
   params.mu = mu;
   params.H = H;
   params.p[0] = p[0];
   params.p[1] = p[1];
   params.v[0] = v[0];
   params.v[1] = v[1];
   params.n = n;
   params.l = l;
   gsl_function f = {&distance_f_st, &params};

   // auxiliary variables
   double htmp;

   // For some reason, bisection method complains that interval [h1,h2] does
   // not straddle 0, so we try to enlarge it a little bit.
   //h1 *= 0.9;
   //h2 /= 0.9;

   if(h2<h1)
   {
	   htmp = h1;
	   h1 = h2;
	   h2 = htmp;
   }

   T = gsl_root_fsolver_brent;
   s = gsl_root_fsolver_alloc (T);
   gsl_root_fsolver_set (s, &f, h1, h2);

   // auxiliary vars
   int status;
   size_t iter = 0;
   double x_lo, x_hi;

   // Find a root of the distance function, i.e. an intersection point of
   // the manifolds (using a bisection method).
   //print_state (iter, s);

    do
      {
	iter++;
	status = gsl_root_fsolver_iterate (s);
	if (status)   /* check if solver is stuck */
	  break;
  
	x_lo = gsl_root_fsolver_x_lower(s);
	x_hi = gsl_root_fsolver_x_upper(s);

	// epsabs=BISECT_TOL, epsrel=0
	status =
	  gsl_root_test_interval (x_lo, x_hi, BISECT_TOL, 0);
        //print_state (iter, s);
      }
    while (status == GSL_CONTINUE && iter < 1000);

    //fprintf (stderr, "status = %s\n", gsl_strerror (status));

    // the root is:
    *h = gsl_root_fsolver_root(s);

    // If bisection did not converge, warn calling function.
    // In this case, the root is updated to the closest zero, 
    // but we don't compute further results, which would be unaccurate (ps,
    // pu, z)
    if(status != GSL_SUCCESS)
       return(2);

   // Compute the following:
   // - point p_s

    p_s[0] = p[0] + (*h) * v[0];
    p_s[1] = p[1] + (*h) * v[1];

   // - intersection point z = P^{-1}(p_s),
   // - t: integration time to reach homoclinic point $z$.

   z[0] = p_s[0];	// x
   z[1] = 0;		// y
   z[2] = p_s[1];	// px
   status=hinv(mu,SEC2,H,z);
   if(status)
   {
      fprintf(stderr, "intersec: error lifting point\n");
      return(1);
   }

   status=prtbp_nl_inv(mu,SEC2,4*n,z,t); 	// $z = P^{-n}(z)$
   if(status)
      {
	 fprintf(stderr, "intersec: error computing intersection point\n");
	 return(1);
      }
   return(0);
}

// name OF FUNCTION: distance_f
// CREDIT: 
//
// DESCRIPTION
// ===========
// We consider the $n$-th iteration under the iterated Poincare map of the
// fundamental segment between the two points $p+h v_u$ and $P(p+h v_u)$.
// Let $p_u = p+h v_u$ be a point in the unstable segment.
//
// We consider the $n$-th iteration of $p_u$ under $\sixmap$, called $q_u$. 
//
// This function computes
//    distance(h) = l - p_x(q_u), 
// where p_x denotes the projection of q_u onto p_x coordinate.
//
// PARAMETERS
// ==========
// h
//    displacement along the unstable direction.
// params
//    mu: mass parameter
//    H: energy value
//    p: fixed point
//    v_u: unstable vector
//    n: number of iterations of the Poincare map.
//    l: axis line
// f
//    On return of the this function, it holds the value
//    distance(h) = l - p_x(q_u).
// 
// RETURN VALUE
// ============
// Returns GSL_SUCCESS to indicate success, or a non-zero error to report an
// error:
//    1: Problems computing the Poincare iterates.
//
// NOTES
// =====

double
distance_f_unst (double h, void *params)
{
   double mu, H;
   double p[2]; 		// fixed point
   double v[2];		// unstable vector

   double n; 			// num. of interations in the unstable dir.
   double l;		// axis line

   double p_u[DIM]; 		// point in the unstable segment

   double t;			// integration time to reach homo. pt.

   // auxiliary vars
   int status;
   double ti;

   mu = ((struct dparams *)params)->mu;
   H = ((struct dparams *)params)->H;

   p[0] = (((struct dparams *)params)->p)[0];
   p[1] = (((struct dparams *)params)->p)[1];

   v[0] = (((struct dparams *)params)->v)[0];
   v[1] = (((struct dparams *)params)->v)[1];

   n = ((struct dparams *)params)->n;
   l = ((struct dparams *)params)->l;

   // Set up point in the unstable segment
   p_u[0] = p[0] + h*v[0];	// x
   p_u[1] = 0;				// y
   p_u[2] = p[1] + h*v[1];	// px
   status=hinv(mu,SEC2,H,p_u);
   if(status)
   {
      fprintf(stderr, "intersec: error lifting point\n");
      return(1);
   }

   // unstable manifold
   status=prtbp_nl(mu,SEC2,4*n,p_u,&t);       // $p_u = P^{n}(p_u)$
   if(status)
   {
      fprintf(stderr, "distance_f: error computing Poincare map\n");
      return(1);
   }

   // Return 
   //    distance(h) = l - p_x(q_u), 
   return l - p_u[2];
}

double
distance_f_st (double h, void *params)
{
   double mu, H;
   double p[2]; 		// fixed point
   double v[2];			// stable vector

   double n; 			// num. of interations in the stable dir.
   double l; 			// axis line

   double p_s[DIM]; 		// point in the stable segment

   double t;			// integration time to reach homo. pt.

   // auxiliary vars
   int status;
   double ti;

   mu = ((struct dparams *)params)->mu;
   H = ((struct dparams *)params)->H;

   p[0] = (((struct dparams *)params)->p)[0];
   p[1] = (((struct dparams *)params)->p)[1];

   v[0] = (((struct dparams *)params)->v)[0];
   v[1] = (((struct dparams *)params)->v)[1];

   n = ((struct dparams *)params)->n;
   l = ((struct dparams *)params)->l;

   // Set up point in the stable segment
   p_s[0] = p[0] + h*v[0];	// x
   p_s[1] = 0;				// y
   p_s[2] = p[1] + h*v[1];	// px
   status=hinv(mu,SEC2,H,p_s);
   if(status)
   {
      fprintf(stderr, "intersec: error lifting point\n");
      return(1);
   }

   // stable manifold
   status=prtbp_nl_inv(mu,SEC2,4*n,p_s,&t);       // $p_s = P^{-n}(p_s)$
   if(status)
   {
      fprintf(stderr, "distance_f: error computing Poincare map\n");
      return(1);
   }

   // Return 
   //    distance(h) = l - p_x(q_u), 
   return l - p_s[2];
}
