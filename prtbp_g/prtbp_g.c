/*! \file
    \brief Poincare map of RTBP in Delaunay coordinates, but integrating the flow in Cartesian.
*/

#include <stdio.h>	    // fprintf
#include <stdlib.h>	    // EXIT_FAILURE
#include <stdbool.h>    // bool
#include <math.h>       // remainder

#include <gsl/gsl_errno.h>	// GSL_SUCCESS
#include <gsl/gsl_roots.h>

#include <frtbp.h>	    // frtbp
#include <rtbp.h>	    // DIM
#include <section.h>	// section_t
#include <cardel.h>	    // cardel

const double TWOPI = 2*M_PI;

const double POINCARE_DEL_CAR_TOL=1.e-16;

/// Integration "step" for prtbp. 
//
// WARNING! This parameter is very delicate. 
// It needs to be appropriately small to get maximum precision for 
// Poincare map. However, setting it too small will increase computation 
// time a lot.
//
// I have used the following values:
// SHORT_TIME_DEL_CAR = 0.001 when computing approximate intersections, 
// and 0.0001 when computing true homoclinic intersections.
const double SHORT_TIME_DEL_CAR=0.001;		

bool onsection_g (section_t sec, double a[DIM]);
bool crossing_fwd_g (section_t sec, double a[DIM], double b[DIM]);
bool crossing_bwd_g (section_t sec, double a[DIM], double b[DIM]);
int inter_g(double mu, section_t sec, double epsabs, 
        double x_del[DIM], double x_car[DIM],
        double t0, double t1, double *t);
double inter_g_f(double t, void *p);

// Parameters to the intersection funtion "inter_g_f"
struct inter_g_f_params;	// Forward declaration

struct inter_g_f_params
{
   double mu; section_t sec;
   double x; double y; double px; double py;
};

/**
 * \remark We do not impose that $x$ is on the section.
 *
 * \remark We do not check that flow at $x$ is transversal to section!!!
 *
 * \remark A point is assumed to be on the Poincare section if it is 
 * within distance POINCARE_DEL_CAR_TOL to the section.
 */

int prtbp_g(double mu, section_t sec, int cuts, double x_del[DIM], 
        double x_car[DIM], double *ti)
{
   double t = 0.0;
   double t_pre;	        /* previous value of time t */
   double x_car_pre[DIM];	/* previous value of point x_car */
   double x_del_pre[DIM];	/* previous value of point x_del */
   int status;
   int i,n;
   double t1;

   n=0;
   while(n!=cuts)
   {
      // Integrate trajectory until it reaches section
      do
      {
         // Save previous value of point "x_del", "x_car" and time "t"
         for(i=0;i<DIM;i++)
         {
            x_del_pre[i]=x_del[i];
            x_car_pre[i]=x_car[i];
         }
         t_pre=t;

         // Integrate for a "short" time t1=SHORT_TIME, short enough so that
         // we can detect crossing of Poincare section.

         // WARNING! Before we used t1=1 as a "short" time, but sometime this
         // was too long...
         status = frtbp(mu,SHORT_TIME_DEL_CAR,x_car);
         t += SHORT_TIME_DEL_CAR;
         if (status != GSL_SUCCESS)
         {
            fprintf(stderr, "prtbp_g: error integrating trajectory\n");
            return(1);
         }
         cardel(x_car,x_del);
      } 
      // while(no crossing of Poincare section)
      while(!(onsection_g(sec,x_del) || 
                  crossing_fwd_g(sec,x_del_pre,x_del))); 
      n++;
      //fprintf(stderr, "DEBUG: befor %d crossing with sect: l=%.15le, g=%.15le\n", n, x_del_pre[0], x_del_pre[2]);
      //fprintf(stderr, "DEBUG: after %d crossing with sect: l=%.15le, g=%.15le\n", n, x_del[0], x_del[2]);
   }
   // point "x_del" is exactly on the section
   // This would be very unlikely...
   if(onsection_g(sec,x_del))
   {
      (*ti)=t;
      return(0);
   }
   // Crossing happened between times t_pre and t. 

   // Restore previous value of point "x_del" and "x_car"
   for(i=0;i<DIM;i++)
   {
      x_del[i]=x_del_pre[i];
      x_car[i]=x_car_pre[i];
   }

   // Intersect trajectory starting at point x with section.
   // WARNING! passing 0 instead of 0.0 gives me trouble?!?!
   if(inter_g(mu, sec, POINCARE_DEL_CAR_TOL, x_del, x_car, 0.0, 
               SHORT_TIME_DEL_CAR, &t1))
   {
          fprintf(stderr, "prtbp_g: error intersecting trajectory with section\n");
          fprintf(stderr, "prtbp_g: giving up...\n");
          /*fprintf(stderr, "x_del_pre: %.15e %.15e %.15e %.15e\n", 
                  x_del_pre[0], x_del_pre[1], x_del_pre[2], x_del_pre[3]);
          fprintf(stderr, "x_car_pre: %.15e %.15e %.15e %.15e\n", 
                  x_car_pre[0], x_car_pre[1], x_car_pre[2], x_car_pre[3]);
          fprintf(stderr, "t: %.15e\n", 
                  SHORT_TIME_DEL_CAR);
          return(1);
                  */
   }
   // Here, point x is on section with tolerance POINCARE_DEL_CAR_TOL. 
   //
   // CAREFUL! Make sure to return a point that is exactly ON the section.
   // If it is slightly before the section, then one more spureous iterate
   // will be counted.
   if(sec==SEC1) x_del[0]=0;
   else if(sec==SEC2) x_del[0]=-M_PI;

   // Set time to reach Poincare section
   (*ti)=t_pre+t1;
   return(0);
}

int prtbp_g_inv(double mu, section_t sec, int cuts, 
        double x_del[DIM], double x_car[DIM], double *ti)
{
   double t = 0.0;
   double t_pre;	        /* previous value of time t */
   double x_car_pre[DIM];	/* previous value of point x_car */
   double x_del_pre[DIM];	/* previous value of point x_del */
   int status;
   int i,n;
   double t1;

   n=0;
   while(n!=cuts)
   {
      // Integrate trajectory backwards until it reaches section
      do
      {
         // Save previous value of point "x_del", "x_car" and time "t"
         for(i=0;i<DIM;i++)
         {
            x_del_pre[i]=x_del[i];
            x_car_pre[i]=x_car[i];
         }
         t_pre=t;

         // Integrate for a "short" time t1=-SHORT_TIME, short enough so that
         // we can detect crossing of Poincare section.

         // WARNING! Before we used t1=1 as a "short" time, but sometime this
         // was too long...
         status = frtbp(mu,-SHORT_TIME_DEL_CAR,x_car);
         t -= SHORT_TIME_DEL_CAR;
         if (status != GSL_SUCCESS)
         {
            fprintf(stderr, "prtbp_g_inv: error integrating trajectory\n");
            return(1);
         }
         cardel(x_car,x_del);
      } 
      // while(no crossing of Poincare section)
      while(!(onsection_g(sec,x_del) || 
                  crossing_fwd_g(sec,x_del_pre,x_del))); 
      n++;
   }
   // point "x_del" is exactly on the section
   // This would be very unlikely...
   if(onsection_g(sec,x_del))
   {
      (*ti)=t;
      return(0);
   }
   // Crossing happened between times t_pre and t. 

   // Restore previous value of point "x_del" and "x_car"
   for(i=0;i<DIM;i++)
   {
      x_del[i]=x_del_pre[i];
      x_car[i]=x_car_pre[i];
   }

   // Intersect trajectory starting at point x with section.
   // Note that (t-t_pre) < 0.
   // WARNING! passing 0 instead of 0.0 gives me trouble?!?!
   if(inter_g(mu, sec, POINCARE_DEL_CAR_TOL, x_del, x_car, 0.0, 
               t-t_pre, &t1))
   {
      fprintf(stderr, "prtbp_g_inv: error intersectig trajectory with section\n");
      return(1);
   }
   // Here, point x is on section with tolerance POINCARE_DEL_CAR_TOL. 

   // Set time to reach Poincare section
   (*ti)=t_pre+t1;
   return(0);
}

/** 
  This function determines if the point A is exactly on the section.

  \param[in] sec 	type of Poincare section (sect = SEC1 or SEC2).
  \param[in] a 		point, 4 coordinates: (l, L, g, G).

  \returns true if point $a$ is exactly on section, false if it is not.
  */
bool onsection_g (section_t sec, double a[DIM])
{
   bool bonsection = false;
   double l;

   l = a[0];

   switch(sec)
   {
      case SEC1 :       // section {l=0}
         {
            bonsection = (fabs(remainder(l,TWOPI))<POINCARE_DEL_CAR_TOL);
            break;
         }
      case SEC2 :       // section {l=\pi}
         {
            bonsection = (fabs(remainder(l-M_PI,TWOPI))<POINCARE_DEL_CAR_TOL);
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
  \param[in] a		First point, 4 coordinates: (l, L, g, G).
  \param[in] b		Second point, 4 coordinates: (l, L, g, G).

  \return 		true if trajectory cuts section, false if it does not.
  */

bool crossing_fwd_g (section_t sec, double a[DIM], double b[DIM])
{
   // auxiliary variables
   double n1,n2;
   bool bcrossing = false;

   switch(sec)
   {
      case SEC1 :       // section {l=0}
         {
             // Since dl/dt>0, it is enough to check if the angle 
             // has been reset from 2\pi to 0.
             // DEBUG bcrossing = (b[0] < a[0]);
             bcrossing = (a[0] < 0 && b[0] > 0);
             break;
         }
      case SEC2 :       // section {l=\pi}
         {
             // Since dl/dt>0, it is enough to check if the angle 
             // has passed from <\pi to >\pi.
             // DEBUG bcrossing = (a[0] < M_PI && b[0] > M_PI);
             bcrossing = (b[0] < a[0]);
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
  \param[in] a		First point, 4 coordinates: (l, L, g, G).
  \param[in] b		Second point, 4 coordinates: (l, L, g, G).

  \return 		true if trajectory cuts section, false if it does not.
  */

bool crossing_bwd_g (section_t sec, double a[DIM], double b[DIM])
{
   // auxiliary variables
   double n1,n2;
   bool bcrossing = false;

   switch(sec)
   {
      case SEC1 :       // section {l=0}
         {
             // Since dl/dt<0, it is enough to check if the angle 
             // has been reset from 0 to 2\pi.
             bcrossing = (b[0] > a[0]);
             break;
         }
      case SEC2 :       // section {l=\pi}
         {
             // Since dl/dt<0, it is enough to check if the angle 
             // has passed from >\pi to <\pi.
             bcrossing = (a[0] > M_PI && b[0] < M_PI);
             break;
         }
   }
   return(bcrossing);
}

double inter_g_f(double t, void *p)
{
   int i, status;
   struct inter_g_f_params *params = (struct inter_g_f_params *)p;
   double mu = (params->mu);
   section_t sec = (params->sec);
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
      fprintf(stderr, "inter_g_f: error computing flow");
      exit(EXIT_FAILURE);
   }

   cardel(pt,pt_del);
   l = pt_del[0];

   switch(sec)
   {
      case SEC1 :       // section {l=0}
         {
            d = remainder(l,TWOPI);
	    //fprintf(stderr,"DEBUG: l=%.15le, d=%.15le\n", l, d);
            break;
         }
      case SEC2 :       // section {l=\pi}
         {
            d = remainder(l-M_PI,TWOPI);
	    //fprintf(stderr,"DEBUG: l=%.15le, d=%.15le\n", l, d);
            break;
         }
   }
   return(d);
}

// name OF FUNCTION: inter_g
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
// "inter_g_f" for the intersection time "t".
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
// x_del
//    Initial point, 4 coordinates: (l, L, g, G). 
//    On return of the function, it holds the intersection point.
// x_car
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
// ERR_MAXITER_DEL_CAR
//    A specified maximum number of iterations of the root finding algorithm
//    has been reached.

const int ERR_MAXITER_DEL_CAR=1;

int inter_g(double mu, section_t sec, double epsabs, 
        double x_del[DIM], double x_car[DIM], 
        double t0, double t1, double *t)
{
    int status;
    int iter = 0, max_iter = 100;
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    double f;
    gsl_function F;
    struct inter_g_f_params params = {mu, sec, x_car[0], x_car[1], x_car[2], x_car[3]};
    *t=0.0;
  
    F.function = &inter_g_f;
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
             "inter_g: solver iteration failed, gsl_errno=%d\n", \
             status);
           return(1);
        }
        *t = gsl_root_fsolver_root (s);
        f=inter_g_f(*t,&params);
        status = gsl_root_test_residual(f, epsabs);
      }
    while (status == GSL_CONTINUE && iter < max_iter);
    gsl_root_fsolver_free (s);

    /*
    if(iter>=max_iter)
    {
       fprintf(stderr, \
               "inter_g: maximum number of iterations reached\n");
       fprintf(stderr, \
               "inter_g: latest residual: %.15e\n",f);
       return(ERR_MAXITER_DEL_CAR);
    }
    */
    // "*t" is the intersection time.
    if(frtbp(mu,*t,x_car))	// compute the intersection point "x"
    {
       fprintf(stderr, "inter_g: error computing intersection point\n");
       exit(EXIT_FAILURE);
    }
    cardel(x_car,x_del);

    if(iter>=max_iter)
    {
       fprintf(stderr, \
               "inter_g: maximum number of iterations reached\n");
       fprintf(stderr, \
               "inter_g: latest residual: %.15e\n",f);
       return(ERR_MAXITER_DEL_CAR);
    }
    return 0;
}
