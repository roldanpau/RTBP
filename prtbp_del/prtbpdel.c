/*! \file
    \brief Poincare map of RTBP in Delaunay coordinates

    $Author: roldan $
    $Date: 2013-03-26 22:23:08 $
*/

// FUNCTIONS
// =========
//
// onsection_del
// -------------
// This function determines if the point x is exactly on the section.
//
// crossing_fwd_del
// ----------------
// This function determines if the trajectory from x to y does cross the
// poincare section sec or not.
//
// crossing_bwd_del
// ----------------
// This function determines if the trajectory from x to y does cross the
// poincare section sec or not.
//
// inter_del
// ---------
// Let S be the Poincare section.
// Given a point "x" whose trajectory intersects the section "S" during the
// time interval (t0,t1), this function computes the intersection time "t"
// such that $flow(t,x)$ is precisely on the section. 
// The intersection point $y = flow(t,x)$ is also returned.
//
// inter_del_f
// -----------
// Function that has to be solved in order to find intersection time.


#include <stdio.h>	// fprintf
#include <stdlib.h>	// EXIT_FAILURE
#include <stdbool.h>	// bool
#include <assert.h>
#include <math.h>	// fmod
#include <gsl/gsl_errno.h>	// GSL_SUCCESS
#include <gsl/gsl_roots.h>
#include <frtbpdel.h>
#include <rtbp.h>	// DIM
#include <section.h>	// section_t, TWOPI

// For $p_0$, I found that asking for 1.e-15 tolerance was too much...
// For iterating $z2_u$ to obtain $z2$, I found that asking for 1.e-14
// tolerance was too much...
const double POINCARE_TOL_DEL=1.e-13;// error bound (tolerance) for Poincare map

bool onsection_del (section_t sec, double x[DIM]);
bool crossing_fwd_del (section_t sec, double x[DIM], double y[DIM]);
bool crossing_bwd_del (section_t sec, double x[DIM], double y[DIM]);
int inter_del(double mu, section_t sec, double epsabs, double x[DIM], double t0, double t1, double *t);
double inter_del_f(double t, void *p);

// Parameters to the intersection funtion "inter_del_f"
struct inter_del_f_params {double mu; section_t sec; double l; double L; double g; double G;};

// NOTES
// =====
// We do not impose that $x$ is on the section.
// We do not check that flow at $x$ is transversal to section!!!
//
// A point is assumed to be on the Poincare section if it is within distance
// POINCARE_TOL_DEL to the section.
//
// Since $l,g$ are angles, we normalize them to \f$[0,2\pi)\f$. 

int prtbp_del(double mu, section_t sec, int cuts, double x[DIM], double *ti)
{
   double t = 0.0;
   double t_pre;	/* previous value of time t */
   double x_pre[DIM];	/* previous value of point x */
   int status;
   int i,n;
   double t1;

   double l,g;
   int q;

   assert(cuts>=0);

   n=0;
   while(n!=cuts)
   {
      // Integrate trajectory until it reaches section
      do
      {
	 // Save previous value of point "x" and time "t"
	 for(i=0;i<DIM;i++)
	    x_pre[i]=x[i];
	 t_pre=t;

	 // Integrate for a "short" time t1=0.1, short enough so that we can
	 // detect crossing of Poincare section.

	 // WARNING! Before we used t1=1 as a "short" time, but sometime this
	 // was too long...
	 status = frtbp_del(mu,0.1,x);
	 t += 0.1;
	 if (status != GSL_SUCCESS)
	 {
	    fprintf(stderr, "prtbp_del: error integrating trajectory\n");
	    return(1);
	 }
      } 
      // while(no crossing of Poincare section)
      while(!(onsection_del(sec,x) || crossing_fwd_del(sec,x,x_pre)));
      n++;
   }
   // point "x" is exactly on the section
   // This would be very unlikely...
   if(onsection_del(sec,x))	
   {
      (*ti)=t;
      return(0);
   }
   // Crossing happened between times t_pre and t. 

   // Restore previous value of point "x"
   for(i=0;i<DIM;i++)
      x[i]=x_pre[i];

   // Intersect trajectory starting at point x with section.
   if(inter_del(mu, sec, POINCARE_TOL_DEL, x, 0, t-t_pre, &t1))
   {
      fprintf(stderr, "prtbp_del: error intersectig trajectory with section\n");
      return(1);
   }
   // Here, point x is on section with tolerance POINCARE_TOL_DEL. 
   // We force x to be exactly on section.
   switch(sec)
   {
      case SEC1 : 	// Poincare section {l=0}
	 {
	    x[0] = 0;		// l
	    break;
	 }
      case SEC2 : 	// Poincare section {l=pi}
	 {
	    x[0] = M_PI;	// l
	    break;
	 }
      case SECg : 	// Poincare section {g=0}
	 {
	    x[2] = 0;	// g
	    break;
	 }
   }

   // Since $l$ is an angle, we normalize it to $[0,2\pi)$.
   l = x[0];
   q = floor(l/TWOPI);
   l -= q*TWOPI;
   x[0]=l;

   // Since $g$ is an angle, we normalize it to $[0,2\pi)$.
   g = x[2];
   q = floor(g/TWOPI);
   g -= q*TWOPI;
   x[2]=g;

   // Set time to reach Poincare section
   (*ti)=t_pre+t1;
   return(0);
}

int prtbp_del_inv(double mu, section_t sec, int cuts, double x[DIM], double *ti)
{
   double t = 0.0;
   double t_pre;	/* previous value of time t */
   double x_pre[DIM];	/* previous value of point x */
   int status;
   int i,n;
   double t1;

   double l,g;
   int q;

   assert(cuts>=0);

   n=0;
   while(n!=cuts)
   {
      // Integrate trajectory until it reaches section
      do
      {
	 // Save previous value of point "x" and time "t"
	 for(i=0;i<DIM;i++)
	    x_pre[i]=x[i];
	 t_pre=t;

	 // Integrate for a "short" time t1=-0.1, short enough so that we can
	 // detect crossing of Poincare section.

	 // WARNING! Before we used t1=-1 as a "short" time, but sometime this
	 // was too long...
	 status = frtbp_del(mu,-0.1,x);
	 t -= 0.1;
	 if (status != GSL_SUCCESS)
	 {
	    fprintf(stderr, "prtbp_del_inv: error integrating trajectory\n");
	    return(1);
	 }
      } 
      // while(no crossing of Poincare section)
      while(!(onsection_del(sec,x) || crossing_bwd_del(sec,x,x_pre)));
      n++;
   }
   // point "x" is exactly on the section
   // This would be very unlikely...
   if(onsection_del(sec,x))	
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
   if(inter_del(mu, sec, POINCARE_TOL_DEL, x, t-t_pre, 0, &t1))
   {
      fprintf(stderr, "prtbp_del_inv: error intersectig trajectory with section\n");
      return(1);
   }
   // Here, point x is on section with tolerance POINCARE_TOL_DEL. 
   // We force x to be exactly on section.
   switch(sec)
   {
      case SEC1 : 	// Poincare section {l=0}
	 {
	    x[0] = 0;		// l
	    break;
	 }
      case SEC2 : 	// Poincare section {l=pi}
	 {
	    x[0] = M_PI;	// l
	    break;
	 }
      case SECg : 	// Poincare section {g=0}
	 {
	    x[2] = 0;	// g
	    break;
	 }
   }

   // Since $l$ is an angle, we normalize it to $[0,2\pi)$.
   l = x[0];
   q = floor(l/TWOPI);
   l -= q*TWOPI;
   x[0]=l;

   // Since $g$ is an angle, we normalize it to $[0,2\pi)$.
   g = x[2];
   q = floor(g/TWOPI);
   g -= q*TWOPI;
   x[2]=g;

   // Set time to reach Poincare section
   (*ti)=t_pre+t1;
   return(0);
}

// name OF FUNCTION: onsection_del
//
// PURPOSE
// =======
// Consider the RTBP in rotating coordinates.
// Let sec be the Poincare section SEC1, SEC2, or SECg corresponding to {l=0},
// {l=\pi}, or {g=0} respectively.
// This function determines if the point x is exactly on the section.
// 
// PARAMETERS
// ==========
// sec
//    type of Poincare section (sect = SEC1, SEC2 or SECg).
// x
//    point, 4 coordinates: (l, L, g, G). 
// 
// RETURN VALUE
// ============
// Returns true if point x is exactly on section, false if it is not.
//
// NOTES
// =====

bool onsection_del (section_t sec, double x[DIM])
{
   bool bonsection = false;

   switch(sec)
   {
      case SEC1 :	// section {l=0}
	 {
         // May be better to say that we are "numerically" on section 
         // if we are within tolerance POINCARE_TOL.
	    bonsection = (fmod(x[0],TWOPI) == 0); 
	    break;
	 }
      case SEC2 :	// section {l=pi}
	 {
	    bonsection = (fmod((x[0]-M_PI),TWOPI) == 0); 
	    break;
	 }
      case SECg :	// section {g=0}
	 {
         // May be better to say that we are "numerically" on section 
         // if we are within tolerance POINCARE_TOL.
	    bonsection = (fmod(x[2],TWOPI) == 0); 
	    break;
	 }
   }
   return(bonsection);
}

// name OF FUNCTION: crossing_fwd_del
//
// PURPOSE
// =======
// Consider the RTBP in rotating coordinates.
// Let sec be the Poincare section SEC1, SEC2 or SECg, corresponding to {l=0}
// {l=\pi} or {g=0} respectively.
// Let x and y be two consecutive points in an orbit (forwards orbit). 
// This function determines if the trajectory from x to y does cross the
// poincare section sec or not.
// 
// PARAMETERS
// ==========
// sec
//    type of Poincare section (sect = SEC1, SEC2 or SECg).
// x
//    First point, 4 coordinates: (l, L, g, G). 
// y
//    Second point, 4 coordinates: (l, L, g, G). 
// 
// RETURN VALUE
// ============
// Returns true if trajectory cuts section, false if it does not.
//
// NOTES
// =====

bool crossing_fwd_del (section_t sec, double x[DIM], double y[DIM])
{
   // auxiliary variables
   double n1,n2;
   bool bCrossing;

   switch(sec)
   {
      case SEC1 :	// section {l=0}
          {
              // Since dl/dt>0, we want to identify [0,2\pi), so we use "floor"
              // function.
              n1 = floor(x[0]/TWOPI);
              n2 = floor(y[0]/TWOPI);
              bCrossing = (n1!=n2);
              break;
          }
      case SEC2 :	// section {l=pi}
          {
              // Since dl/dt>0, we want to identify [0,2\pi), so we use "floor"
              // function.
              n1 = floor((x[0]-M_PI)/TWOPI);
              n2 = floor((y[0]-M_PI)/TWOPI);
              bCrossing = (n1!=n2);
              break;
          }
      case SECg :	// section {g=0}
          {
              // Since dg/dt<0, we want to identify (-2\pi,0], so we use "ceil"
              // function.
              n1 = ceil(x[2]/TWOPI);
              n2 = ceil(y[2]/TWOPI);
              // Inevitably, g will jump by almost TWOPI when (x,y) changes
              // from the 2nd quadrant to the 3rd (see cardel.c).
              // We need to exclude this case as a "fake" crossing of section.
              bCrossing = ((n1!=n2) && !(fabs(x[2]-y[2])>(TWOPI-0.2)));
              break;
          }
   }
   return(bCrossing);
}

// name OF FUNCTION: crossing_bwd_del
//
// PURPOSE
// =======
// Consider the RTBP in rotating coordinates.
// Let sec be the Poincare section SEC1, SEC2 or SECg, corresponding to {l=0}
// {l=\pi} or {g=0} respectively.
// Let x and y be two consecutive points in an orbit (backwards orbit). 
// This function determines if the trajectory from x to y does cross the
// poincare section sec or not.
// 
// PARAMETERS
// ==========
// sec
//    type of Poincare section (sect = SEC1, SEC2 or SECg).
// x
//    First point, 4 coordinates: (l, L, g, G). 
// y
//    Second point, 4 coordinates: (l, L, g, G). 
// 
// RETURN VALUE
// ============
// Returns true if trajectory cuts section, false if it does not.
//
// NOTES
// =====

bool crossing_bwd_del (section_t sec, double x[DIM], double y[DIM])
{
   // auxiliary variables
   double n1,n2;
   bool bCrossing;

   switch(sec)
   {
      case SEC1 :	// section {l=0}
          {
              // Since dl/dt<0, we want to identify (-2\pi,0), so we use "ceil"
              // function.
              n1 = ceil(x[0]/TWOPI);
              n2 = ceil(y[0]/TWOPI);
              bCrossing = (n1!=n2);
              break;
          }
      case SEC2 :	// section {l=pi}
          {
              // Since dl/dt<0, we want to identify (-2\pi,0), so we use "ceil"
              // function.
              n1 = ceil((x[0]-M_PI)/TWOPI);
              n2 = ceil((y[0]-M_PI)/TWOPI);
              bCrossing = (n1!=n2);
              break;
          }
      case SECg :	// section {g=0}
          {
              // Since dg/dt>0, we want to identify [0,2\pi), so we use "floor"
              // function.
              n1 = floor(x[2]/TWOPI);
              n2 = floor(y[2]/TWOPI);
              // Inevitably, g will jump by almost TWOPI when (x,y) changes
              // from the 2nd quadrant to the 3rd (see cardel.c).
              // We need to exclude this case as a "fake" crossing of section.
              bCrossing = ((n1!=n2) && !(fabs(x[2]-y[2])>(TWOPI-0.1)));
              break;
          }
   }
   return(bCrossing);
}

// name OF FUNCTION: inter_del
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
// "inter_del_f" for the intersection time "t".
//
// PARAMETERS
// ==========
// mu
//    mass parameter of the RTBP
// sec
//    type of Poincare section = {SEC1,SEC2,SECg}
// epsabs
//    maximum desired error bound (tolerance) for intersection.
//    A point $p=(l,L,g,G)$ is considered to intersect the section if
//       dist(p,sec) < epsabs.
// x
//    Initial point, 4 coordinates: (l, L, g, G). 
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
// ERR_MAXITER_DEL
//    A specified maximum number of iterations of the root finding algorithm
//    has been reached.
//
// NOTES
// =====
//
// CALLS TO: inter_del_f, frtbp

const int ERR_MAXITER_DEL=1;

int inter_del(double mu, section_t sec, double epsabs, double x[DIM], double t0, double t1, double *t)
{
    int status;
    int iter = 0, max_iter = 100;
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    double f;
    gsl_function F;
    struct inter_del_f_params params = {mu, sec, x[0], x[1], x[2], x[3]};
    *t=0.0;
  
    F.function = &inter_del_f;
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
		 "inter_del: solver iteration failed, gsl_errno=%d\n", status);
	   return(1);
	}
	*t = gsl_root_fsolver_root (s);
	f=inter_del_f(*t,&params);
	status = gsl_root_test_residual(f, epsabs);
      }
    while (status == GSL_CONTINUE && iter < max_iter);
    gsl_root_fsolver_free (s);

    if(iter>=max_iter)
    {
       fprintf(stderr, "inter_del: maximum number of iterations reached\n");
       return(ERR_MAXITER_DEL);
    }
    // "*t" is the intersection time.
    if(frtbp_del(mu,*t,x))	// compute the intersection point "x"
    {
       fprintf(stderr, "inter_del: error computing intersection point\n");
       exit(EXIT_FAILURE);
    }
    return 0;
}

double inter_del_f(double t, void *p)
{
   int i, status;
   struct inter_del_f_params *params = (struct inter_del_f_params *)p;
   double mu = (params->mu);
   section_t sec = (params->sec);
   double l = (params->l); 
   double L = (params->L); 
   double g = (params->g);
   double G = (params->G);
   double pt[DIM];

   // auxiliary variables
   double d;	// ditance to section

   pt[0]=l; 
   pt[1]=L; 
   pt[2]=g;
   pt[3]=G;

   status=frtbp_del(mu,t,pt);	// flow(t,pt)
   if(status!=0)
   {
      fprintf(stderr, "inter_del_f: error computing flow");
      exit(EXIT_FAILURE);
   }

   switch(sec)
   {
      case SEC1 :
	 {
	    d = remainder(pt[0],TWOPI);
	    break;
	 }
      case SEC2 : 
	 {
	    d = remainder((pt[0]-M_PI),TWOPI);
	    break;
	 }
      case SECg :
	 {
	    d = remainder(pt[2],TWOPI);
	    break;
	 }
   }
   return(d);
}
