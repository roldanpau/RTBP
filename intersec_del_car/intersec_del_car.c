/*! \file intersec_del_car.c
    \brief Intersection of Invariant Manifolds

    $Author: roldan $
    $Date: 2013-03-11 11:27:15 $
*/

#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_roots.h>

#include <math.h>   // remainder

#include <rtbp.h>		// DIM
#include <hinv.h>
#include <cardel.h>
#include <prtbp_2d.h>		// prtbp_2d_inv
#include <prtbp_del_car.h>
#include <lift.h>
#include <utils_module.h>	// dblcpy

/// Tolerance (precision) for bisection method.
//
// Note that this is the requested precision for the interval 
// (h_1, h_2) on the local unstable manifold, NOT for the 
// homoclinic point itself!
const double BISECT_TOL=1.e-14;

// Parameters to distance_f_unst and distance_f_st functions.
struct dparams
{
   double mu;		// mass parameter
   section_t sec;   // Poincare section 
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

 fprintf (stderr, "iter = %3u, [%.7f, %.7f], root = %.7f, err(est) = %.15f\n",
	 (unsigned int)iter,
	 x_lo, 
	 x_hi,
	 r,
	 x_hi-x_lo);
}

/**
  Take a point in the local unstable manifold and iterate it n times by 
  the Poincare map.

  Consider the 2D map 
  \f$\mathcal{P}: S \to S\f$, where S is the Poincare section SEC1 or SEC2, 
  corresponding to \f$\{l=0\}\f$ or \f$\{l=\pi\}\f$.

  Let $p$ be a hyperbolic fixed point for \f$\mathcal{P}\f$.
  For definiteness, we assume that $p$ is located above the $p_x=0$ axis.

  Assume that \f$\lambda\f$ is the unstable eigenvalue, with \f$\lambda>1\f$.
  Let $v$ be the unstable eigenvector for the eigenvalue \f$\lambda\f$. 
  For definiteness, we assume that $v=(x,p_x)$ points "to the right", i.e. we
  assume that the first component of $v$ is $x>0$, but we could use the other
  branch of the manifold.
  Let $W^u(p)$ be the unstable manifold of $p$.

  This function computes $z$, the $n$-th iteration under the iterated 
  Poincare map of the point $p_u = p+h_u v$ located in the unstable segment.

  \param[in] mu         mass parameter for the RTBP
  \param[in] sec        Poincare section: sec={SEC1,SEC2}
  \param[in] H          energy value

  \param[in] n
  number of desired iterations by the Poincare map in the unstable direction.

  \param[in] p_u[2]
  Point in the unstable segment (Cartesian coords)

  \param[out] t
  On exit, it contains integration time to reach z from p_u.
  
  \param[out] z_del[DIM]	
  On exit, it contains homoclinic point z = P^n(p_u) in Delaunay.

  \param[out] z_car[DIM]	
  On exit, it contains homoclinic point z = P^n(p_u) in Cartesian.

  \param[out] z_u[DIM]
  Point in local unstable manifold of the pendulum (Delaunay coords). This will
  be needed in outer_circ.
  
  \param[out] z_u_car[DIM]
  Point in local unstable manifold of the pendulum (Cartesian coords). This will
  be needed in outer_circ.
  
  \returns 
  a non-zero error code to indicate an error and 0 to indicate
  success.

  \retval 1 	Problems computing the Poincare iterates.
*/

// NOTES
// =====
// For the moment, we work with the 3:1 resonant family of periodic orbits.

int iterate_del_car_unst(double mu, section_t sec, double H, int n, 
        double p_u[2], double *t, double z_del[DIM], double z_car[DIM],
	double z_u[DIM], double z_u_car[DIM])
{
    // auxiliary variables
    int status;

   // Lift point from \R^2 to \R^4
   status=lift(mu,SEC2,H,1,p_u,z_car);
   if(status)
   {
      perror("intersec_del_car: error lifting point");
      return(1);
   }

   // Transform point to Delaunay coordinates
   cardel(z_car,z_del);

   // Flow the point to Delaunay section before the iteration.
   status=prtbp_del_car(mu,sec,1,z_del,z_car,t);  // $q_u = P^{n}(p_u)$
   if(status)
   {
      fprintf(stderr, "intersec_del_car: error flowing point to Delaunay sec\n");
      return(1);
   }

   // Return z_u
   dblcpy(z_u,z_del, DIM);
   dblcpy(z_u_car,z_car, DIM);

   // unstable manifold
   status=prtbp_del_car(mu,sec,n,z_del,z_car,t);  // $q_u = P^{n}(p_u)$
   if(status)
   {
      fprintf(stderr, "intersec_del_car: error computing Poincare map\n");
      //return(1);
   }
   return(0);
}

/**
  Intersection of unstable invariant manifold with symmetry line.

  Consider the 2D map 
  \f$\mathcal{P}: S \to S\f$, where S is the Poincare section SEC1 or SEC2, 
  corresponding to \f$\{l=0\}\f$ or \f$\{l=\pi\}\f$, which is
  assumed to be reversible with respect to the symmetry line $g=0$ and $g=\pi$.

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
  \param[in] sec        Poincare section: sec={SEC1,SEC2}
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
  
  \param[out] z_del[DIM]	
  On exit, it contains homoclinic point z = P^n(p_u) in Delaunay.

  \param[out] z_car[DIM]	
  On exit, it contains homoclinic point z = P^n(p_u) in Cartesian.

  \param[out] z_u[DIM]
  Point in local unstable manifold of the appropriate pendulum (Delaunay
  coords). This will be needed in outer_circ.
  
  \param[out] z_u_car[DIM]
  Point in local unstable manifold of the appropriate pendulum (Cartesian
  coords). This will be needed in outer_circ.
  
  \returns 
  a non-zero error code to indicate an error and 0 to indicate
  success.

  \retval 1 	Problems computing the Poincare iterates.
  \retval 2 	Bisection procedure did not converge
*/

// NOTES
// =====
// For the moment, we work with the 3:1 resonant family of periodic orbits.

int intersec_del_car_unst(double mu, section_t sec, double H, double p[2], 
        double v[2], double lambda, int n, double h1, double h2, double l,
        double *h, double p_u[2], double *t, double z_del[DIM], 
        double z_car[DIM], double z_u[DIM], double z_u_car[DIM])
{
   const gsl_root_fsolver_type *T;
   gsl_root_fsolver *s;

   struct dparams params;
   params.mu = mu;
   params.sec = sec;
   params.H = H;
   params.p[0] = p[0];
   params.p[1] = p[1];
   params.v[0] = v[0];
   params.v[1] = v[1];
   params.n = n;
   params.l = l;
   gsl_function f = {&distance_f_unst, &params};

   // Bisection method sometimes complains that interval [h1,h2] does not
   // straddle 0. The reason is that approxint_del_car sets h1 and h2 from
   // points in the linear segment l, whereas intersec_del_car uses h1 and h2
   // to recover those endpoints of the segment u_i, so there is a small
   // discrepancy.
   // Thus we enlarge [h1,h2] a little bit to account for this discrepancy.
   // NOTE: The number 0.9 is crucial! We tried 0.999 and did not work...
   h1 *= 0.97;
   h2 /= 0.97;

   T = gsl_root_fsolver_brent;
   s = gsl_root_fsolver_alloc (T);

   // For first branch, h1<h2, so...
   gsl_root_fsolver_set (s, &f, h1, h2);

   // For second branch, h2<h1, so...
   //gsl_root_fsolver_set (s, &f, h2, h1);

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

  // - Point in local unstable manifold of the appropriate pendulum (Delaunay
  // coords). This will be needed in outer_circ.
  // - Point in local unstable manifold of the appropriate pendulum (Cartesian
  // coords). This will be needed in outer_circ.
   // - intersection point z = P(p_u),
   // - t: integration time to reach homoclinic point $z$.
   status=iterate_del_car_unst(mu,sec,H,n,p_u,t,z_del,z_car,z_u,z_u_car);
   if(status)
   {
      perror("intersec_del_car: error iterating point");
      return(1);
   }
   return(0);
}

/**
  Intersection of stable invariant manifold with symmetry line.

  Exactly as \ref intersec_del_car_unst.
  */

int intersec_del_car_st(double mu, double H, double p[2], double v[2], 
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

   // For some reason, bisection method complains that interval [h1,h2] does
   // not straddle 0, so we try to enlarge it a little bit.
   h1 *= 0.9;
   h2 /= 0.9;

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

   z[0] = p_s[0];
   z[1] = p_s[1];
   status=prtbp_2d_inv(mu,SEC2,H,2*n,z,t); 	// $z = P^{-n}(z)$
   if(status)
      {
	 fprintf(stderr, "intersec_del_car: error computing intersection point\n");
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
   section_t sec;   // Poincare section
   double p[2]; 		// fixed point
   double v[2];		// unstable vector

   double n; 			// num. of interations in the unstable dir.
   double l;		// axis line

   double p_u[2]; 		// point in the unstable segment

   double t;			// integration time to reach homo. pt.

   double d;            // distance(h) = l - p_x(q_u).

   // auxiliary vars
   int status;
   double z_car[DIM]; 		// homoclinic point (cartesian)
   double z_del[DIM]; 		// homoclinic point (Delaunay)
   double z_u[DIM]; 		// point in the local unstable manifold (Delaunay)
   double z_u_car[DIM]; 	// point in the local unstable manifold (Cartesian)

   mu = ((struct dparams *)params)->mu;
   sec = ((struct dparams *)params)->sec;
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

   status=iterate_del_car_unst(mu,sec,H,n,p_u,&t,z_del,z_car,z_u,z_u_car);
   if(status)
   {
      perror("intersec_del_car: error iterating point");
      return(1);
   }

   if(sec==SEC1 || sec==SEC2)
   {
       // Return distance(h) = l - g(q_u), 
       d=remainder(z_del[2]-l,2*M_PI);
   }
   else if(sec==SECg || sec==SECg2)
   {
       // Return distance(h) = l - l(q_u), 
       d=remainder(z_del[0]-l,2*M_PI);
       //fprintf(stderr, "distance=%.16e\n", d);
   }
   return d;
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
   status=prtbp_2d_inv(mu,SEC2,H,2*n,p_s,&t);       // $p_s = P^{-n}(p_s)$
   if(status)
   {
      fprintf(stderr, "distance_f: error computing Poincare map\n");
      return(1);
   }

   // Return 
   //    distance(h) = l - p_x(q_u), 
   return l - p_s[1];
}
