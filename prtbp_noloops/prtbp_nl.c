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
const double SHORT_TIME_NL=0.001;		///< integration "step" for prtbp_nl

int inter_nl(double mu, section_t sec, double epsabs, double x[DIM], double t0,
      double t1, double *t);
double inter_nl_f(double t, void *p);

// Parameters to the intersection funtion "inter_nl_f"
struct inter_nl_f_params;	// Forward declaration

struct inter_nl_f_params
{
   double mu; section_t sec;
   double x; double y; double px; double py;
};

/** 
  This function determines if the point A is exactly on the section.

  \param[in] sec 	type of Poincare section (sect = SEC1 or SEC2).
  \param[in] a 		point, 4 coordinates: (x, y, p_x, p_y).

  \returns true if point $a$ is exactly on section, false if it is not.
  */
bool onsection_nl (section_t sec, double a[DIM])
{
   bool bonsection = false;

   double x = a[0];
   double y = a[1];
   double py = a[3];
   double vy = py-x;

   /*
   switch(sec)
   {
      case SEC1 :       // section {y=0, v_y>0}
         {
            bonsection = (y == 0 && vy > 0);
            break;
         }
      case SEC2 :       // section {y=0, v_y<0}
         {
            bonsection = (y == 0 && vy < 0);
            break;
         }
   }
   */
   bonsection = (y == 0);
   return(bonsection);
}

/**
  Let $a$ and $b$ be two consecutive points in an orbit.
  This function determines if the trajectory from $a$ to $b$ does cross the
  poincare section sec or not.

  \param[in] sec 	type of Poincare section (sec = SEC1 or SEC2).
  \param[in] a		First point, 4 coordinates: (x, y, p_x, p_y).
  \param[in] b		Second point, 4 coordinates: (x, y, p_x, p_y).

  \return 		true if trajectory cuts section, false if it does not.
  */
// NOTES
// =====
// In fact, we should check if $v_y>0$ (or $v_y<0$) at the point $c$ that's
// precisely on the section, not at the point $b$.
// However, we don't expect $v_y$ to change much between these two points.

bool crossing_nl (section_t sec, double a[DIM], double b[DIM])
{
   double x = a[0];
   double py = a[3];
   double vy = py-x;

   // auxiliary variables
   double n1,n2;

   bool bcrossing = false;

   /*
   switch(sec)
   {
      case SEC1 :       // section {y=0, v_y>0}
         {
            bcrossing = (a[1]*b[1]<0 && vy>0);
            break;
         }
      case SEC2 :       // section {y=0, v_y<0}
         {
            bcrossing = (a[1]*b[1]<0 && vy<0);
            break;
         }
   }
   */
   bcrossing = (a[1]*b[1]<0);
   return(bcrossing);
}

int onemorecut_nl_1(double mu, section_t sec, double x[DIM], double *t, 
		double x_pre[DIM], double *t_pre)
{
	int i;
	int status;

    int sign = ((x[0]<0) ? -1:1);

	do {
	 
		// Save previous value of point "x" and time "t"
		for(i=0;i<DIM;i++)
			x_pre[i]=x[i];
		*(t_pre)=(*t);
	 
		// Integrate for a "short" time t1=SHORT_TIME_NL, short enough so that
		// we can detect crossing of Poincare section.

		// WARNING! Before we used t1=1 as a "short" time, but sometime this
		// was too long...
		status = frtbp(mu,SHORT_TIME_NL,x);
		(*t) += SHORT_TIME_NL;
		if (status != GSL_SUCCESS)
		{
			fprintf(stderr, "prtbp_nl: error integrating trajectory\n");
			return(1);
		}
	} 
	// while(no crossing of Poincare section)
	while(!(onsection_nl(sec,x) || crossing_nl(sec,x_pre,x)) || (sign*x[0]>0)); 
}

int onemorecut_nl_2(double mu, section_t sec, double x[DIM], double *t, 
		double x_pre[DIM], double *t_pre)
{
	int i;
	int status;

	do {
	 
		// Save previous value of point "x" and time "t"
		for(i=0;i<DIM;i++)
			x_pre[i]=x[i];
		*(t_pre)=(*t);
	 
		// Integrate for a "short" time t1=SHORT_TIME_NL, short enough so that
		// we can detect crossing of Poincare section.

		// WARNING! Before we used t1=1 as a "short" time, but sometime this
		// was too long...
		status = frtbp(mu,SHORT_TIME_NL,x);
		(*t) += SHORT_TIME_NL;
		if (status != GSL_SUCCESS)
		{
			fprintf(stderr, "prtbp_nl: error integrating trajectory\n");
			return(1);
		}
	} 
	// while(no crossing of Poincare section)
	while(!(onsection_nl(sec,x) || crossing_nl(sec,x_pre,x))); 
}

// NOTES
// =====
// We do not impose that $x$ is on the section.
// We do not check that flow at $x$ is transversal to section!!!
//
// A point is assumed to be on the Poincare section if it is within distance
// POINCARE_TOL_NL to the section.
//

int prtbp_nl(double mu, section_t sec, int cuts, double x[DIM], double *ti)
{
   double t = 0.0;
   double t_pre;	/* previous value of time t */
   double x_pre[DIM];	/* previous value of point x */
   int i,n;
   double t1;

   // auxiliary variables
   double x1[DIM], x2[DIM];
   double tfirst, tsecond;
   double t_pre1, t_pre2;
   double x_pre1[DIM], x_pre2[DIM];

   tfirst = 0.0;
   n=0;
   while(n!=cuts)
   {
	   tfirst=t;
	   dblcpy(x1,x,DIM);

	  // Integrate trajectory until it reaches section
	  onemorecut_nl_1(mu,sec,x1,&tfirst,x_pre1,&t_pre1);
	  tsecond=tfirst;
	  dblcpy(x2,x1,DIM);
	  onemorecut_nl_2(mu,sec,x2,&tsecond,x_pre2,&t_pre2);

	  if(x_pre1[0]*x_pre2[0]<0)
	  {
		  t=tfirst;
		  dblcpy(x,x1,DIM);
		  dblcpy(x_pre,x_pre1,DIM);
		  t_pre=t_pre1;
	  }
	  if(x_pre1[0]*x_pre2[0]>0)
	  {
          t=tsecond;
          dblcpy(x,x2,DIM);
          dblcpy(x_pre,x_pre2,DIM);
          t_pre=t_pre2;
	  }
      n++;
   }
   // point "x" is exactly on the section
   // This would be very unlikely...
   if(onsection_nl(sec,x))
   {
      (*ti)=t;
      return(0);
   }
   // Crossing happened between times t_pre and t. 

   // Restore previous value of point "x"
   for(i=0;i<DIM;i++)
      x[i]=x_pre[i];

   // Intersect trajectory starting at point x with section.
   // WARNING! passing 0 instead of 0.0 gives me trouble?!?!
   if(inter_nl(mu, sec, POINCARE_TOL_NL, x, 0.0, t-t_pre, &t1))
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
// We do not check that flow at $x$ is transversal to section!!!
//
// A point is assumed to be on the Poincare section if it is within distance
// POINCARE_TOL_NL to the section.
//
// On successful return of this function, the point $x$ is exactly on the
// section, i.e. we set coordinate $y$ exactly equal to zero.

int prtbp_nl_inv(double mu, section_t sec, int cuts, double x[DIM], double *ti)
{
   double t = 0.0;
   double t_pre;	/* previous value of time t */
   double x_pre[DIM];	/* previous value of point x */
   int sign;		/* sign of previous intersection with x axis */
   int status;
   int i,n;
   double t1;

   // Save sign of previous intersection with x axis
   sign = (x[0]>0 ? +1 : -1);

   n=0;
   while(n!=cuts)
   {
      // Integrate trajectory backwards until it reaches section
      do
      {
	 // Save previous value of point "x" and time "t"
	 for(i=0;i<DIM;i++)
	    x_pre[i]=x[i];
	 t_pre=t;

	 // Integrate for a "short" time t1=-SHORT_TIME_NL, short enough so that
	 // we can detect crossing of Poincare section.

	 // WARNING! Before we used t1=1 as a "short" time, but sometime this
	 // was too long...
	 status = frtbp(mu,-SHORT_TIME_NL,x);
	 t -= SHORT_TIME_NL;
	 if (status != GSL_SUCCESS)
	 {
	    fprintf(stderr, "prtbp_inv: error integrating trajectory\n");
	    return(1);
	 }
      } 
      // while(no crossing of Poincare section)
      while(!(onsection_nl(sec,x) || crossing_nl(sec,x_pre,x))); 
      n++;
   }
   // point "x" is exactly on the section
   // This would be very unlikely...
   if(onsection_nl(sec,x))
   {
      (*ti)=t;
      return(0);
   }
   // Crossing happened between times t_pre and t. 

   // Restore previous value of point "x"
   for(i=0;i<DIM;i++)
      x[i]=x_pre[i];

   // Intersect trajectory starting at point x with section.
   // Note that (t-t_pre) < 0.
   // WARNING! passing 0 instead of 0.0 gives me trouble?!?!
   if(inter_nl(mu, sec, POINCARE_TOL_NL, x, t-t_pre, 0.0, &t1))
   {
      fprintf(stderr, "prtbp_inv: error intersectig trajectory with section\n");
      return(1);
   }
   // Here, point x is on section with tolerance POINCARE_TOL_NL_DEL. 
   // We force x to be exactly on section.
   x[1] = 0;    // y

   // Set time to reach Poincare section
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

int inter_nl(double mu, section_t sec, double epsabs, double x[DIM], double t0,
      double t1, double *t)
{
    int status;
    int iter = 0, max_iter = 100;
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    double f;
    gsl_function F;
    struct inter_nl_f_params params = {mu, sec, x[0], x[1], x[2], x[3]};
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
       fprintf(stderr, "inter_nl: maximum number of iterations reached\n");
       return(ERR_MAXITER_NL);
    }
    // "*t" is the intersection time.
    if(frtbp(mu,*t,x))	// compute the intersection point "x"
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
   //section_t sec = (params->sec);	// not used
   double x = (params->x); 
   double y = (params->y); 
   double px = (params->px);
   double py = (params->py);
   double pt[DIM];

   pt[0]=x; 
   pt[1]=y; 
   pt[2]=px;
   pt[3]=py;

   status=frtbp(mu,t,pt);	// flow(t,pt)
   if(status!=0)
   {
      fprintf(stderr, "inter_nl_f: error computing flow");
      exit(EXIT_FAILURE);
   }
   return(pt[1]);
}
