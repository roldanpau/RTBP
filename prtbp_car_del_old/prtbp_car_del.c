/*! \file
    \brief Poincare map of RTBP in Delaunay coordinates, but integrating the flow in Cartesian.
*/

#include <stdio.h>	// fprintf
#include <stdlib.h>	// EXIT_FAILURE
#include <stdbool.h>    // bool
#include <math.h>    // remainder

#include <gsl/gsl_errno.h>	// GSL_SUCCESS
#include <gsl/gsl_roots.h>

#include <frtbp.h>	// frtbp
#include <rtbp.h>	// DIM
#include <section.h>	// section_t
#include <cardel.h>	// cardel

const double TWOPI = 2*M_PI;

const double POINCARE_CAR_DEL_TOL=1.e-16;
const double SHORT_TIME=1;		///< integration "step" for prtbp

int inter_car_del(double mu, section_t sec, double epsabs, double x[DIM], double t0,
      double t1, double *t);
double inter_car_del_f(double t, void *p);

// Parameters to the intersection funtion "inter_car_del_f"
struct inter_car_del_f_params;	// Forward declaration

struct inter_car_del_f_params
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
bool onsection_car_del (section_t sec, double a[DIM])
{
   bool bonsection = false;
   double b[DIM];   // (l,L,g,G)
   double l;

   cardel(a,b);
   l = b[0];

   switch(sec)
   {
      case SEC1 :       // section {l=0}
         {
            bonsection = (fabs(remainder(l,TWOPI))<POINCARE_CAR_DEL_TOL);
            break;
         }
      case SEC2 :       // section {l=\pi}
         {
            bonsection = (fabs(remainder(l-M_PI,TWOPI))<POINCARE_CAR_DEL_TOL);
            break;
         }
   }
   return(bonsection);
}

/**
  Let $a$ and $b$ be two consecutive points in an orbit (forwards orbit).
  This function determines if the trajectory from $a$ to $b$ does cross the
  poincare section sec or not.

  \param[in] sec 	type of Poincare section (sec = SEC1 or SEC2).
  \param[in] a		First point, 4 coordinates: (x, y, p_x, p_y).
  \param[in] b		Second point, 4 coordinates: (x, y, p_x, p_y).

  \return 		true if trajectory cuts section, false if it does not.
  */

bool crossing_fwd_car_del (section_t sec, double a[DIM], double b[DIM])
{
   // auxiliary variables
   double n1,n2;
   double x[DIM], y[DIM];

   cardel(a,x);
   cardel(b,y);

   bool bcrossing = false;

   switch(sec)
   {
      case SEC1 :       // section {l=0}
         {
             // Since dl/dt>0, it is enough to check if the angle 
             // has been reset from 2\pi to 0.
             bcrossing = (y[0] < x[0]);
             break;
         }
      case SEC2 :       // section {l=\pi}
         {
             // Since dl/dt>0, it is enough to check if the angle 
             // has passed from <\pi to >\pi.
             bcrossing = (x[0] < M_PI && y[0] > M_PI);
             break;
         }
   }
   return(bcrossing);
}

/**
  Let $a$ and $b$ be two consecutive points in an orbit (backwards orbit).
  This function determines if the trajectory from $a$ to $b$ does cross the
  poincare section sec or not.

  \param[in] sec 	type of Poincare section (sec = SEC1 or SEC2).
  \param[in] a		First point, 4 coordinates: (x, y, p_x, p_y).
  \param[in] b		Second point, 4 coordinates: (x, y, p_x, p_y).

  \return 		true if trajectory cuts section, false if it does not.
  */

bool crossing_bwd_car_del (section_t sec, double a[DIM], double b[DIM])
{
   // auxiliary variables
   double n1,n2;
   double x[DIM], y[DIM];

   cardel(a,x);
   cardel(b,y);

   bool bcrossing = false;

   switch(sec)
   {
      case SEC1 :       // section {l=0}
         {
             // Since dl/dt<0, it is enough to check if the angle 
             // has been reset from 0 to 2\pi.
             bcrossing = (y[0] > x[0]);
             break;
         }
      case SEC2 :       // section {l=\pi}
         {
             // Since dl/dt<0, it is enough to check if the angle 
             // has passed from >\pi to <\pi.
             bcrossing = (x[0] > M_PI && y[0] < M_PI);
             break;
         }
   }
   return(bcrossing);
}

// NOTES
// =====
// We do not impose that $x$ is on the section.
// We do not check that flow at $x$ is transversal to section!!!
//
// A point is assumed to be on the Poincare section if it is within distance
// POINCARE_CAR_DEL_TOL to the section.

int prtbp_car_del(double mu, section_t sec, int cuts, double x[DIM], 
        double *ti)
{
   double t = 0.0;
   double t_pre;	/* previous value of time t */
   double x_pre[DIM];	/* previous value of point x */
   int status;
   int i,n;
   double t1;

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

         // Integrate for a "short" time t1=SHORT_TIME, short enough so that
         // we can detect crossing of Poincare section.

         // WARNING! Before we used t1=1 as a "short" time, but sometime this
         // was too long...
         status = frtbp(mu,SHORT_TIME,x);
         t += SHORT_TIME;
         if (status != GSL_SUCCESS)
         {
            fprintf(stderr, "prtbp: error integrating trajectory\n");
            return(1);
         }
      } 
      // while(no crossing of Poincare section)
      while(!(onsection_car_del(sec,x) || 
                  crossing_fwd_car_del(sec,x_pre,x))); 
      n++;
   }
   // point "x" is exactly on the section
   // This would be very unlikely...
   if(onsection_car_del(sec,x))
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
   if(inter_car_del(mu, sec, POINCARE_CAR_DEL_TOL, x, 0.0, t-t_pre, &t1))
   {
      fprintf(stderr, "prtbp: error intersectig trajectory with section\n");
      return(1);
   }
   // Here, point x is on section with tolerance POINCARE_CAR_DEL_TOL. 

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
// POINCARE_CAR_DEL_TOL to the section.

int prtbp_inv(double mu, section_t sec, int cuts, double x[DIM], double *ti)
{
   double t = 0.0;
   double t_pre;	/* previous value of time t */
   double x_pre[DIM];	/* previous value of point x */
   int status;
   int i,n;
   double t1;

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

         // Integrate for a "short" time t1=-SHORT_TIME, short enough so that
         // we can detect crossing of Poincare section.

         // WARNING! Before we used t1=1 as a "short" time, but sometime this
         // was too long...
         status = frtbp(mu,-SHORT_TIME,x);
         t -= SHORT_TIME;
         if (status != GSL_SUCCESS)
         {
            fprintf(stderr, "prtbp_inv: error integrating trajectory\n");
            return(1);
         }
      } 
      // while(no crossing of Poincare section)
      while(!(onsection_car_del(sec,x) || 
                  crossing_bwd_car_del(sec,x_pre,x))); 
      n++;
   }
   // point "x" is exactly on the section
   // This would be very unlikely...
   if(onsection_car_del(sec,x))
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
   if(inter_car_del(mu, sec, POINCARE_CAR_DEL_TOL, x, t-t_pre, 0.0, &t1))
   {
      fprintf(stderr, "prtbp_inv: error intersectig trajectory with section\n");
      return(1);
   }
   // Here, point x is on section with tolerance POINCARE_CAR_DEL_TOL. 

   // Set time to reach Poincare section
   (*ti)=t_pre+t1;
   return(0);
}

// name OF FUNCTION: inter_car_del
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
// "inter_car_del_f" for the intersection time "t".
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
//    dist(p,sec) < epsabs.
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
// ERR_MAXITER_CAR_DEL
//    A specified maximum number of iterations of the root finding algorithm
//    has been reached.

const int ERR_MAXITER_CAR_DEL=1;

int inter_car_del(double mu, section_t sec, double epsabs, double x[DIM], 
        double t0, double t1, double *t)
{
    int status;
    int iter = 0, max_iter = 100;
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    double f;
    gsl_function F;
    struct inter_car_del_f_params params = {mu, sec, x[0], x[1], x[2], x[3]};
    *t=0.0;
  
    F.function = &inter_car_del_f;
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
             "inter_car_del: solver iteration failed, gsl_errno=%d\n", \
             status);
           return(1);
        }
        *t = gsl_root_fsolver_root (s);
        f=inter_car_del_f(*t,&params);
        status = gsl_root_test_residual(f, epsabs);
      }
    while (status == GSL_CONTINUE && iter < max_iter);
    gsl_root_fsolver_free (s);

    if(iter>=max_iter)
    {
       fprintf(stderr, \
               "inter_car_del: maximum number of iterations reached\n");
       return(ERR_MAXITER_CAR_DEL);
    }
    // "*t" is the intersection time.
    if(frtbp(mu,*t,x))	// compute the intersection point "x"
    {
       fprintf(stderr, "inter: error computing intersection point\n");
       exit(EXIT_FAILURE);
    }
    return 0;
}

double inter_car_del_f(double t, void *p)
{
   int i, status;
   struct inter_car_del_f_params *params = (struct inter_car_del_f_params *)p;
   double mu = (params->mu);
   section_t sec = (params->sec);	// not used
   double x = (params->x); 
   double y = (params->y); 
   double px = (params->px);
   double py = (params->py);
   double pt[DIM];

   // auxiliary variables
   double pt_del[DIM];  // point in Delaunay coordinates
   double l;
   double d;            // distance to section

   pt[0]=x; 
   pt[1]=y; 
   pt[2]=px;
   pt[3]=py;

   status=frtbp(mu,t,pt);	// flow(t,pt)
   if(status!=0)
   {
      fprintf(stderr, "inter_car_del_f: error computing flow");
      exit(EXIT_FAILURE);
   }

   cardel(pt,pt_del);
   l = pt_del[0];

   switch(sec)
   {
      case SEC1 :       // section {l=0}
         {
            d = remainder(l,TWOPI);
            break;
         }
      case SEC2 :       // section {l=\pi}
         {
            d = remainder(l-M_PI,TWOPI);
            break;
         }
   }
   return(d);
}
