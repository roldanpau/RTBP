// ===================================
// Intersection of Invariant Manifolds
// ===================================
// FILE:          $RCSfile: intersecdel.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2012-06-06 11:07:40 $
//
// FUNCTIONS
// =========
//
// intersec_del_unst, intersec_del_st
// ----------------------------------
// Consider the Restricted Three Body Problem in Delaunay variables. 
// Assume that energy is fixed: $H(l,L,g,G)=\bar H$.
// Let SEC2= \{l=\pi\}$ be the Poincare section, and $\sixmap$ be the
// ``iterated'' 2D Poincare map.
// Let $p$ be a fixed point of $\sixmap$ associated to the 3:1 resonant
// periodic orbit of the flow, and let $W^u(p)$ be the associated unstable
// manifold.
// This function computes the ``primary'' intersection point $z=(l,G)$ with
// the axis line defined by
//    \[ g = l. \]
// Due to the symmetry of the manifolds in this section, it is enough to look
// for the first intersection of the unstable manifold $W^u(p)$ with the
// reversibility line (the $G$ axis of the poincare section).

#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_roots.h>

#include <rtbp.h>		// DIM
#include <prtbpdel_2d.h>	// prtbp_del_2d, prtbp_del_2d_inv

// Tolerance (precision) for bisection method 
const double BISECT_TOL=1.e-15;

struct dparams
{
   double mu;		// mass parameter
   double H;		// energy value
   double p[2];		// fixed point
   double v[2];		// unstable/stable vector
   double n;		// num. of iterations in the unstable/stable dir.
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

 fprintf (stderr, "iter = %3u, [%.7f, %.7f], root = %.7f, err(est) = %.15f\n",
	 iter,
	 x_lo, 
	 x_hi,
	 r,
	 x_hi-x_lo);
}

// name OF FUNCTION: intersec_del_unst
//
// PURPOSE
// =======
// Consider the Restricted Three Body Problem in Delaunay variables. 
// Assume that energy is fixed: $H(l,L,g,G)=\bar H$.
// Let SEC2= \{l=\pi\}$ be the Poincare section, and $\sixmap$ be the
// ``iterated'' 2D Poincare map.
// Let $p$ be a fixed point of $\sixmap$ associated to the 3:1 resonant
// periodic orbit of the flow, and let $W^u(p)$ be the associated unstable
// manifold.
// This function computes the ``primary'' intersection point $z=(l,G)$ with
// the axis line defined by
//    \[ g = l. \]
// Due to the symmetry of the manifolds in this section, it is enough to look
// for the first intersection of the unstable manifold $W^u(p)$ with the
// reversibility line (the $G$ axis of the poincare section).
//
// We consider the $n$-th iteration under the iterated Poincare map of the
// fundamental segment between the two points $p+h v_u$ and $P(p+h v_u)$.
// Let $p_u = p+h_u v_u$ be a point in the unstable segment.
//
// We are provided (by the program approxint) with an interval $(h_1,h_2)$
// that is guaranteed to contain the root $h^*$ such that 
//    \[ g(P(p+h^* v_u))=l. \]
//
// Then we use a standard numerical method (bisection-like 1-dimensional root
// finding) to find $h^*$.
//
// The function that we solve (find a zero of) is
//    distance(h_u) = l - g(P(p+h_u v_u)).
//
// PARAMETERS
// ==========
// mu
//    mass parameter for the RTBP
// H
//    energy value
// p[2]
//    fixed point p=(g,G)
// v_u[2]
//    unstable vector
// lambda_u
//    unstable eigenvalue
// n
//    number of desired iterations by the Poincare map in the unstable
//    direction.
//    This is furnished by program "approxint".
// h_1, h_2
//    Small increment in the direction of v_u (initial interval for bisection
//    method).
//    This is furnished by program "approxint".
// l
//    Axis line
// h
//    On exit, it contains the true root found by numerical method.
// p_u[2]
//    On exit, it contains point in the unstable segment such that
//       g(P(p_u)) = 0.
// tu
//    On exit, it contains integration time to reach z from p_u.
// z[2]	
//    On exit, it contains homoclinic point z = P^n(p_u).
// 
// RETURN VALUE
// ============
// Returns a non-zero error code to indicate an error and 0 to indicate
// success:
//    1: Problems computing the Poincare iterates.
//    2: Bisection procedure did not converge
//
// NOTES
// =====
// For the moment, we work with the 3:1 resonant family of periodic orbits.
//
// We assume that the manifolds lie on the SEC2 section {l=pi}.
//
// We ask for precision of BISECT_TOL in the homoclinc point.

int intersec_del_unst(double mu, double H, double p[2], double v[2], 
      double lambda, int n, double h1, double h2, double l,
      double *h, double p_u[2], double *t, double z[2])
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

   T = gsl_root_fsolver_brent;
   s = gsl_root_fsolver_alloc (T);
   gsl_root_fsolver_set (s, &f, h1, h2);

   // auxiliary vars
   int status;
   size_t iter = 0;
   double x_lo, x_hi;

   // Find a root of the distance function, i.e. an intersection point of
   // the manifolds (using a bisection method).
   print_state (iter, s);

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
        print_state (iter, s);
      }
    while (status == GSL_CONTINUE && iter < 1000);

    fprintf (stderr, "status = %s\n", gsl_strerror (status));

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

   z[0] = p_u[0];
   z[1] = p_u[1];
   status=prtbp_del_2d(mu,SEC2,H,3*n,z,t); 	// $z = P^n(z)$
   if(status)
      {
	 fprintf(stderr, "intersec_del: error computing intersection point\n");
	 return(1);
      }
   return(0);
}

int intersec_del_st(double mu, double H, double p[2], double v[2], 
      double lambda, int n, double h1, double h2, double l,
      double *h, double p_s[2], double *t, double z[2])
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

   T = gsl_root_fsolver_brent;
   s = gsl_root_fsolver_alloc (T);
   gsl_root_fsolver_set (s, &f, h1, h2);

   // auxiliary vars
   int status;
   size_t iter = 0;
   double x_lo, x_hi;

   // Find a root of the distance function, i.e. an intersection point of
   // the manifolds (using a bisection method).
   print_state (iter, s);

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
        print_state (iter, s);
      }
    while (status == GSL_CONTINUE && iter < 1000);

    fprintf (stderr, "status = %s\n", gsl_strerror (status));

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

   z[0] = p_s[0];
   z[1] = p_s[1];
   status=prtbp_del_2d_inv(mu,SEC2,H,3*n,z,t); 	// $z = P^{-n}(z)$
   if(status)
      {
	 fprintf(stderr, "intersec_del: error computing intersection point\n");
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
// Let $p_u = p+h_u v_u$ be a point in the unstable segment.
//
// We consider the $n$-th iteration of $p_u$ under $\sixmap$, called $q_u$. 
//
// This function computes
//    distance(h_u) = g(q_u), 
// i.e. the projection of q_u onto g coordinate.
//
// PARAMETERS
// ==========
// h_u
//    displacement along the unstable direction.
// params
//    mu: mass parameter
//    H: energy value
//    p: fixed point
//    v_u: unstable vector
//    n: number of iterations of the Poincare map.
// f
//    On return of the this function, it holds the value
//    distance(h_u) = g(q_u).
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

   double p_u[2]; 		// point in the unstable segment

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
   p_u[0] = p[0] + h*v[0];
   p_u[1] = p[1] + h*v[1];

   // unstable manifold
   status=prtbp_del_2d(mu,SEC2,H,3*n,p_u,&t);       // $p_u = P^{n}(p_u)$
   if(status)
   {
      fprintf(stderr, "distance_f: error computing Poincare map\n");
      return(1);
   }

   // Return 
   //    distance(h) = l - g(q_u), 
   return l - p_u[0];
}

double
distance_f_st (double h, void *params)
{
   double mu, H;
   double p[2]; 		// fixed point
   double v[2];			// stable vector

   double n; 			// num. of interations in the stable dir.
   double l; 			// axis line

   double p_s[2]; 		// point in the stable segment

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
   p_s[0] = p[0] + h*v[0];
   p_s[1] = p[1] + h*v[1];

   // stable manifold
   status=prtbp_del_2d_inv(mu,SEC2,H,3*n,p_s,&t);       // $p_s = P^{-n}(p_s)$
   if(status)
   {
      fprintf(stderr, "distance_f: error computing Poincare map\n");
      return(1);
   }

   // Return 
   //    distance(h) = l - g(q_u), 
   return l - p_s[0];
}
