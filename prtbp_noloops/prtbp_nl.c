/*! \file
    \brief Poincare map of the Restricted Three Body Problem

    $Author: roldan $
    $Date: 2013-03-26 22:22:14 $
*/

#include <stdio.h>	// fprintf
#include <stdlib.h>	// EXIT_FAILURE
#include <stdbool.h>    // bool

#include <gsl/gsl_errno.h>	// GSL_SUCCESS
#include <gsl/gsl_roots.h>

#include <frtbp.h>	// frtbp
#include <rtbp.h>	// DIM

#include <section.h>	// section_t
#include <utils_module.h>	// dblcpy

const double POINCARE_TOL_NL=1.e-16;
const double TANGENT_TOL_NL=1.e-6;     ///< tolerance for tangent condition
const double SHORT_TIME_NL=0.01;		///< integration "step" for prtbp_nl

int inter_nl(double mu, double epsabs, double x[DIM], double t0, double t1,
		bool fwd, double *t);
double inter_nl_f(double t, void *p);

// Parameters to the intersection funtion "inter_nl_f"
struct inter_nl_f_params;	// Forward declaration

struct inter_nl_f_params
{
   double mu; double x; double y; double px; double py; bool fwd;
};

/** 
  This function determines if the flow is tangent to the Poincare 
  section at the point A.

  \param[in] sec 	type of Poincare section (sect = SEC1 or SEC2).
  \param[in] a 		point, 4 coordinates: (x, y, p_x, p_y).

  \returns true the flow is tangent, false if it is not.
  */
bool tangent_nl(section_t sec, double a[DIM])
{
   double x = a[0];
   double y = a[1];
   double py = a[3];
   double vy = py-x;

   // This works for both sections:
   if(fabs(y)<TANGENT_TOL_NL && fabs(vy)<TANGENT_TOL_NL)
   {
       fprintf(stderr, "tranverse_nl: Flow is tangent to section!!\n");
       return(true);
   }
   return(false);
}

/** 
  This function determines if the point A is exactly on the section.

  \param[in] sec 	type of Poincare section (sect = SEC1 or SEC2).
  \param[in] a 		point, 4 coordinates: (x, y, p_x, p_y).

  \param[in,out] sign
    (pointer to) sign of previous intersection with x axis (sign = +1 if x>0
    or sign = -1 if x<0).

  \returns true if point $a$ is exactly on section, false if it is not.
  */
bool onsection_nl (section_t sec, double a[DIM], int *sign)
{
   bool bonsection = false;

   double x = a[0];
   double y = a[1];
   double py = a[3];
   double vy = py-x;

   // TWO CONSECUTIVE ITERATES must NOT lie both to the right or to the left
   // of the origin, so they must have different sign.
   bonsection = ((y == 0) && (*sign)*x<0);

   // If we intersected the x axis, determine if we crossed to the left
   // or to the right of the origin.
   if(y == 0)
      *sign = (x>0 ? +1 : -1);

   return(bonsection);
}

/**
  Let $a$ and $b$ be two consecutive points in an orbit.
  This function determines if the trajectory from $a$ to $b$ does cross the
  poincare section sec or not.

  \param[in] sec 	type of Poincare section (sec = SEC1 or SEC2).
  \param[in] a		First point, 4 coordinates: (x, y, p_x, p_y).
  \param[in] b		Second point, 4 coordinates: (x, y, p_x, p_y).

  \param[in,out] sign
    (pointer to) sign of previous intersection with x axis (sign = +1 if x>0
    or sign = -1 if x<0).

  \return 		true if trajectory cuts section, false if it does not.
  */
// NOTES
// =====
// In fact, we should check if $v_y>0$ (or $v_y<0$) at the point $c$ that's
// precisely on the section, not at the point $b$.
// However, we don't expect $v_y$ to change much between these two points.

bool crossing_nl (section_t sec, double a[DIM], double b[DIM], int *sign)
{
   double x = a[0];
   double py = a[3];
   double vy = py-x;

   bool bcrossing = false;

   // TWO CONSECUTIVE ITERATES must NOT lie both to the right or to the left
   // of the origin, so they must have different sign.
   bcrossing = ((a[1]*b[1]<0) && (*sign)*x<0);

   // If we have intersected the x axis, determine if we crossed to the left
   // or to the right of the origin.
   if(a[1]*b[1]<0)
      *sign = (x>0 ? +1 : -1);

   return(bcrossing);
}

// NOTES
// =====
// We do not impose that $x$ is on the section.
//
// A point is assumed to be on the Poincare section if it is within distance
// POINCARE_TOL_NL to the section.
//
// Possible improvement: In frtbp, check if the flow becomes tangential to
// section at some point.  If it does, return a flag to prtbp, which sould act
// accordingly.

int prtbp_nl(double mu, section_t sec, int cuts, double x[DIM], double *ti)
{
   double t = 0.0;
   double t_pre;	/* previous value of time t */
   double t_pre2;	/* previous value of time t */
   double x_pre[DIM];	/* previous value of point x */
   double x_pre2[DIM];	/* previous value of point x */
   double x2[DIM];	
   double t2;
   int sign_pre;		/* sign of previous intersection with x axis */
   int sign_cur;		/* sign of current intersection with x axis */
   int sign_nxt;		/* sign of next intersection with x axis */
   int status;
   int n;
   double t1;

   if(tangent_nl(sec,x))
   {
       perror("Flow is tangent to section. Cannot compute Poincare map!\n");
       exit(EXIT_FAILURE);
   }

   n=0;
   while(n<cuts)
   {
	   // Save sign of previous intersection with x axis
	   sign_pre = (x[0]>0 ? +1 : -1);

	  // Integrate trajectory until it crosses x axis
	  do
	  {
		 // Save previous value of point "x" and time "t"
		 dblcpy(x_pre, x, DIM);
		 t_pre=t;

		 // Integrate for a "short" time t1=SHORT_TIME_NL, short enough so that
		 // we can detect crossing of Poincare section.

		 // WARNING! Before we used t1=1 as a "short" time, but sometime this
		 // was too long...
		 status = frtbp(mu,SHORT_TIME_NL,x);
		 t += SHORT_TIME_NL;
		 if (status != GSL_SUCCESS)
		 {
			fprintf(stderr, "prtbp_nl: error integrating trajectory\n");
			return(1);
		 }
	  } 
	  // while(no crossing of x axis)
	  while(!(x[1] == 0 || x_pre[1]*x[1] < 0)); 

	   // Save sign of current intersection with x axis
	   sign_cur = (x[0]>0 ? +1 : -1);

	  dblcpy(x2, x, DIM);
	  t2 = t;
	  // Integrate trajectory until it crosses x axis one more time
	  do
	  {
		 dblcpy(x_pre2, x2, DIM);
		 t_pre2=t2;

		 // Integrate for a "short" time t1=SHORT_TIME_NL, short enough so that
		 // we can detect crossing of Poincare section.

		 // WARNING! Before we used t1=1 as a "short" time, but sometime this
		 // was too long...
		 status = frtbp(mu,SHORT_TIME_NL,x2);
		 t2 += SHORT_TIME_NL;
		 if (status != GSL_SUCCESS)
		 {
			fprintf(stderr, "prtbp_nl: error integrating trajectory\n");
			return(1);
		 }
	  } 
	  // while(no crossing of x axis)
	  while(!(x2[1] == 0 || x_pre2[1]*x2[1] < 0)); 

	   // Save sign of next intersection with x axis
	   sign_nxt = (x2[0]>0 ? +1 : -1);

	  if((sign_pre!=sign_cur && sign_cur!=sign_nxt) || 
			  (sign_pre==sign_cur && sign_cur==sign_nxt)) n++;
	  //else
	//	  fprintf(stderr, "prtbp_nl: skipping cut with x axis...\n");
   }

   // point "x" is exactly on the section
   // This would be very unlikely...
   if(x[1] == 0)
   {
      (*ti)=t;
      return(0);
   }
   // Crossing happened between times t_pre and t. 

   // Restore previous value of point "x"
   dblcpy(x, x_pre, DIM);

   // Intersect trajectory starting at point x with section.
   // Integrate fwd.
   // WARNING! passing 0 instead of 0.0 gives me trouble?!?!
   if(inter_nl(mu, POINCARE_TOL_NL, x, 0.0, t-t_pre, true, &t1))
   {
      fprintf(stderr, "prtbp_nl: error intersectig trajectory with section\n");
      return(1);
   }
   // Here, point x is on section with tolerance POINCARE_TOL_NL_DEL. 
   // We force x to be exactly on section.
   x[1] = 0;    // y

   // Set time to reach Poincare section
   (*ti)=t_pre+t1;
   return(0);
}

// NOTES
// =====
// We do not impose that $x$ is on the section.
//
// A point is assumed to be on the Poincare section if it is within distance
// POINCARE_TOL_NL to the section.
//
// Possible improvement: In frtbp, check if the flow becomes tangential to
// section at some point.  If it does, return a flag to prtbp, which sould act
// accordingly.

int prtbp_nl_inv(double mu, section_t sec, int cuts, double x[DIM], double *ti)
{
   double t = 0.0;
   double t_pre;	/* previous value of time t */
   double t_pre2;	/* previous value of time t */
   double x_pre[DIM];	/* previous value of point x */
   double x_pre2[DIM];	/* previous value of point x */
   double x2[DIM];	
   double t2;
   int sign_pre;		/* sign of previous intersection with x axis */
   int sign_cur;		/* sign of current intersection with x axis */
   int sign_nxt;		/* sign of next intersection with x axis */
   int status;
   int n;
   double t1;

   if(tangent_nl(sec,x))
   {
       perror("Flow is tangent to section. Cannot compute Poincare map!\n");
       exit(EXIT_FAILURE);
   }

   n=0;
   while(n<cuts)
   {
	   // Save sign of previous intersection with x axis
	   sign_pre = (x[0]>0 ? +1 : -1);

	  // Integrate trajectory until it crosses x axis
	  do
	  {
		 // Save previous value of point "x" and time "t"
		 dblcpy(x_pre, x, DIM);
		 t_pre=t;

		 // Integrate for a "short" time t1=-SHORT_TIME_NL, short enough so that
		 // we can detect crossing of Poincare section.

		 // WARNING! Before we used t1=1 as a "short" time, but sometime this
		 // was too long...
		 status = frtbp(mu,-SHORT_TIME_NL,x);
		 t -= SHORT_TIME_NL;
		 if (status != GSL_SUCCESS)
		 {
			fprintf(stderr, "prtbp_nl_inv: error integrating trajectory\n");
			return(1);
		 }
	  } 
	  // while(no crossing of x axis)
	  while(!(x[1] == 0 || x_pre[1]*x[1] < 0)); 

	   // Save sign of current intersection with x axis
	   sign_cur = (x[0]>0 ? +1 : -1);

	  dblcpy(x2, x, DIM);
	  t2 = t;
	  // Integrate trajectory until it crosses x axis one more time
	  do
	  {
		 dblcpy(x_pre2, x2, DIM);
		 t_pre2=t2;

		 // Integrate for a "short" time t1=-SHORT_TIME_NL, short enough so that
		 // we can detect crossing of Poincare section.

		 // WARNING! Before we used t1=1 as a "short" time, but sometime this
		 // was too long...
		 status = frtbp(mu,-SHORT_TIME_NL,x2);
		 t2 -= SHORT_TIME_NL;
		 if (status != GSL_SUCCESS)
		 {
			fprintf(stderr, "prtbp_nl_inv: error integrating trajectory\n");
			return(1);
		 }
	  } 
	  // while(no crossing of x axis)
	  while(!(x2[1] == 0 || x_pre2[1]*x2[1] < 0)); 

	   // Save sign of next intersection with x axis
	   sign_nxt = (x2[0]>0 ? +1 : -1);

	  if((sign_pre!=sign_cur && sign_cur!=sign_nxt) || 
			  (sign_pre==sign_cur && sign_cur==sign_nxt)) n++;
	  //else
	//	  fprintf(stderr, "prtbp_nl: skipping cut with x axis...\n");
   }

   // point "x" is exactly on the section
   // This would be very unlikely...
   if(x[1] == 0)
   {
      (*ti)=t;
      return(0);
   }
   // Crossing happened between times t_pre and t. 

   // Restore previous value of point "x"
   dblcpy(x, x_pre, DIM);

   // Intersect trajectory starting at point x with section.
   // Integrate bwd.
   // Note that (t-t_pre) < 0.
   // WARNING! passing 0 instead of 0.0 gives me trouble?!?!
   if(inter_nl(mu, POINCARE_TOL_NL, x, 0.0, -(t-t_pre), false, &t1))
   {
      fprintf(stderr, "prtbp_nl_inv: error intersectig trajectory with section\n");
      return(1);
   }
   // Here, point x is on section with tolerance POINCARE_TOL_NL_DEL. 
   // We force x to be exactly on section.
   x[1] = 0;    // y

   // Set time to reach Poincare section
   t1 = -t1;
   (*ti)=t_pre+t1;
   return(0);
}

// name OF FUNCTION: inter_nl
// CREDIT: 
//
// PURPOSE
// =======
// Let "sec" be the given Poincare section.
// Given a point "x" whose trajectory intersects the section "sec" during the
// time interval (t0,t1), this function computes the intersection time "t"
// such that $flow(t,x)$ is precisely on the section. 
// The intersection point $y = flow(t,x)$ is also returned.
//
// We use a (1-dimensional) root finding algorithm to solve the function
// "inter_nl_f" for the intersection time "t".
//
// PARAMETERS
// ==========
// mu
//    mass parameter of the RTBP
// sec
//    type of Poincare section = {SEC1,SEC2}
// epsabs
//    maximum desired error bound (tolerance) for intersection.
//    A point $p=(x,y,p_x,p_y)$ is considered to intersect the section if
//       |y| < epsabs.
// x
//    Initial point, 4 coordinates: (X, Y, P_X, P_Y). 
//    On return of the function, it holds the intersection point.
// t0,t1
//    Lower and upper bounds for the time of intersection. 
//    The caller must guarantee that solution lies in between (t0,t1).
// t
//    Pointer to intersection time.
//    On return of the this function, *t holds the intersection time.
// 
// RETURN VALUE
// ============
// Returns a non-zero error code to indicate an error and 0 to indicate
// success.
// If an error is encountered, the function returns a non-zero value:
//
// ERR_MAXITER_NL
//    A specified maximum number of iterations of the root finding algorithm
//    has been reached.
//
// NOTES
// =====
//
// CALLS TO: inter_nl_f, frtbp

const int ERR_MAXITER_NL=1;

int inter_nl(double mu, double epsabs, double x[DIM], double t0,
      double t1, bool fwd, double *t)
{
    int status;
    int iter = 0, max_iter = 100;
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    double f;
    gsl_function F;
    struct inter_nl_f_params params = {mu, x[0], x[1], x[2], x[3], fwd};
    *t=0.0;
  
    F.function = &inter_nl_f;
    F.params = &params;
  
    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc (T);
    gsl_root_fsolver_set (s, &F, t0, t1);

    do
      {
	iter++;
	status = gsl_root_fsolver_iterate (s);
	if (status)
	{
	   fprintf(stderr, \
		 "inter_nl: solver iteration failed, gsl_errno=%d\n", status);
	   return(1);
	}
	*t = gsl_root_fsolver_root (s);
	f=inter_nl_f(*t,&params);
	status = gsl_root_test_residual(f, epsabs);
      }
    while (status == GSL_CONTINUE && iter < max_iter);
    gsl_root_fsolver_free (s);

    if(iter>=max_iter)
    {
       fprintf(stderr, "inter_nl: maximum number of iterations reached.\n");
       fprintf(stderr, "inter_nl: Last residual: %e\n", f);
       return(ERR_MAXITER_NL);
    }
    // "*t" is the intersection time.
	// Compute the intersection point "x"
	if(fwd)
		status = frtbp(mu,*t,x);
	else
		status = frtbp(mu,-(*t),x);
    if(status)	
    {
       fprintf(stderr, "inter_nl: error computing intersection point\n");
       exit(EXIT_FAILURE);
    }
    return 0;
}

double inter_nl_f(double t, void *p)
{
   int i, status;
   struct inter_nl_f_params *params = (struct inter_nl_f_params *)p;
   double mu = (params->mu);
   double x = (params->x); 
   double y = (params->y); 
   double px = (params->px);
   double py = (params->py);
   bool fwd = (params->fwd);
   double pt[DIM];

   pt[0]=x; 
   pt[1]=y; 
   pt[2]=px;
   pt[3]=py;

	// flow(t,pt)
	if(fwd)
		status = frtbp(mu,t,pt);
	else
		status = frtbp(mu,-t,pt);
   if(status!=0)
   {
      fprintf(stderr, "inter_nl_f: error computing flow");
      exit(EXIT_FAILURE);
   }
   return(pt[1]);
}
