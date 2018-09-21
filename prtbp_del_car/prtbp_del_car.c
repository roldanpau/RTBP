/*! \file
    \brief Poincare map of RTBP in Delaunay coordinates, but integrating the
    flow in Cartesian.
*/

#include <stdio.h>	    // fprintf
#include <stdlib.h>	    // EXIT_FAILURE
#include <stdbool.h>    // bool
#include <math.h>       // remainder

#include <gsl/gsl_errno.h>	// GSL_SUCCESS
#include <gsl/gsl_roots.h>

#include <frtbp.h>	    // frtbp
#include <rtbp.h>	    // DIM
#include <section.h>	// TWOPI, section_t
#include <cardel.h>	    // cardel
#include <prtbpdel.h>	// crossing_fwd_del, crossing_bwd_del

#include <utils_module.h>	    // dblcpy

// Unfortunately, this code does not get as much precision as prtbp_del. 
// Empirically, we get precisions close to 1.e-8. 
// The reason is that the problem of evaluating the function $g(phi(t))$ is
// highly sensitive: small changes to the argument $t$ lead to large changes in
// the solution $g(phi(t))$.
const double POINCARE_DEL_CAR_TOL=1.e-16;

/// Integration "step" for prtbp. 
//
// WARNING! This parameter is very delicate. 
// It needs to be appropriately small to get maximum precision for 
// Poincare map. However, setting it too small will increase computation 
// time a lot.
//
// I have used the following vals:
// SHORT_TIME_DEL_CAR = 0.001 when computing approximate intersections, 
// and 0.0001 when computing true homoclinic intersections.

// PRG (04/04/2018): const double SHORT_TIME_DEL_CAR=0.001;		
const double SHORT_TIME_DEL_CAR=0.1;

bool onsection_del_car (section_t sec, double a[DIM]);
int inter_del_car(double mu, section_t sec, double epsabs, 
        double x_del[DIM], double x_car[DIM],
        double t0, double t1, double *t);
double inter_del_car_f(double t, void *p);

// Parameters to the intersection funtion "inter_del_car_f"
struct inter_del_car_f_params;	// Forward declaration

struct inter_del_car_f_params
{
   double mu; section_t sec;
   double x; double y; double px; double py;
};

double WrapPosNegPI(double fAng)
{
	return fmod(fAng + M_PI, TWOPI) - M_PI;
}

/**
 * \remark We do not impose that $x$ is on the section.
 *
 * \remark We do not check that flow at $x$ is transversal to section!!!
 *
 * \remark A point is assumed to be on the Poincare section if it is 
 * within distance POINCARE_DEL_CAR_TOL to the section.
 */

int prtbp_del_car(double mu, section_t sec, int cuts, double x_del[DIM], 
        double x_car[DIM], double *ti)
{
   double t = 0.0;
   double t_pre;	        /* previous value of time t */
   double x_car_pre[DIM];	/* previous value of point x_car */
   double x_del_pre[DIM];	/* previous value of point x_del */
   int status;
   int i,n;
   double t1;
   double g;

   /* Normalize g between (-pi,pi]
   g=x_del[2];
   if(g<=-M_PI || g>M_PI)
   {
	   x_del[2] = WrapPosNegPI(g);
   } */

   n=0;
   while(n!=cuts)
   {
      // Integrate trajectory until it reaches section
      do
      {
         // Save previous value of point "x_del", "x_car" and time "t"
         dblcpy(x_del_pre, x_del, DIM);
         dblcpy(x_car_pre, x_car, DIM);
         t_pre=t;

         // Integrate for a "short" time t1=SHORT_TIME, short enough so that
         // we can detect crossing of Poincare section.

         // WARNING! Before we used t1=1 as a "short" time, but sometime this
         // was too long...
         status = frtbp(mu,SHORT_TIME_DEL_CAR,x_car);
         t += SHORT_TIME_DEL_CAR;
         if (status != GSL_SUCCESS)
         {
            fprintf(stderr, "prtbp_del_car: error integrating trajectory\n");
            return(1);
         }
         cardel(x_car,x_del);
      } 
      // while(no crossing of Poincare section)
      while(!(onsection_del_car(sec,x_del) || 
                  crossing_fwd_del(sec,x_del_pre,x_del))); 
      n++;
      //fprintf(stderr, "DEBUG: befor %d crossing with sect: l=%.15le, g=%.15le\n", n, x_del_pre[0], x_del_pre[2]);
      //fprintf(stderr, "DEBUG: after %d crossing with sect: l=%.15le, g=%.15le\n", n, x_del[0], x_del[2]);
   }
   // point "x_del" is exactly on the section
   // This would be very unlikely...
   if(onsection_del_car(sec,x_del))
   {
      (*ti)=t;
      return(0);
   }
   // Crossing happened between times t_pre and t. 

   // Restore previous value of point "x_del" and "x_car"
   dblcpy(x_del, x_del_pre, DIM);
   dblcpy(x_car, x_car_pre, DIM);

   // Intersect trajectory starting at point x with section.
   // WARNING! passing 0 instead of 0.0 gives me trouble?!?!
   if(inter_del_car(mu, sec, POINCARE_DEL_CAR_TOL, x_del, x_car, 0.0, 
               t-t_pre, &t1))
   {
          fprintf(stderr, "prtbp_del_car: error intersecting trajectory with section\n");
          fprintf(stderr, "prtbp_del_car: giving up...\n");
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
   else if(sec==SEC2) x_del[0]=-M_PI;   // Since dl/dt>0
   else if(sec==SECg) x_del[2]=0;
   else if(sec==SECg2) x_del[2]=M_PI;   // Since dg/dt<0

   // Set time to reach Poincare section
   (*ti)=t_pre+t1;
   return(0);
}

int prtbp_del_car_inv(double mu, section_t sec, int cuts, 
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
         dblcpy(x_del_pre, x_del, DIM);
         dblcpy(x_car_pre, x_car, DIM);
         t_pre=t;

         // Integrate for a "short" time t1=-SHORT_TIME, short enough so that
         // we can detect crossing of Poincare section.

         // WARNING! Before we used t1=1 as a "short" time, but sometime this
         // was too long...
         status = frtbp(mu,-SHORT_TIME_DEL_CAR,x_car);
         t -= SHORT_TIME_DEL_CAR;
         if (status != GSL_SUCCESS)
         {
            fprintf(stderr, "prtbp_del_car_inv: error integrating trajectory\n");
            return(1);
         }
         cardel(x_car,x_del);
      } 
      // while(no crossing of Poincare section)
      while(!(onsection_del_car(sec,x_del) || 
                  crossing_bwd_del(sec,x_del_pre,x_del))); 
      n++;
   }
   // point "x_del" is exactly on the section
   // This would be very unlikely...
   if(onsection_del_car(sec,x_del))
   {
      (*ti)=t;
      return(0);
   }
   // Crossing happened between times t_pre and t. 

   // Restore previous value of point "x_del" and "x_car"
   dblcpy(x_del, x_del_pre, DIM);
   dblcpy(x_car, x_car_pre, DIM);

   // Intersect trajectory starting at point x with section.
   // Note that (t-t_pre) < 0.
   // WARNING! passing 0 instead of 0.0 gives me trouble?!?!
   if(inter_del_car(mu, sec, POINCARE_DEL_CAR_TOL, x_del, x_car, t-t_pre, 0.0,
               &t1))
   {
          fprintf(stderr, "prtbp_del_car: error intersecting trajectory with section\n");
          fprintf(stderr, "prtbp_del_car: giving up...\n");
   }
   // Here, point x is on section with tolerance POINCARE_DEL_CAR_TOL. 
   //
   // CAREFUL! Make sure to return a point that is exactly ON the section.
   // If it is slightly before the section, then one more spureous iterate
   // will be counted.
   if(sec==SEC1) x_del[0]=0;
   else if(sec==SEC2) x_del[0]=M_PI;   // Since dl/dt>0
   else if(sec==SECg) x_del[2]=0;
   else if(sec==SECg2) x_del[2]=-M_PI;   // Since dg/dt<0

   // Set time to reach Poincare section
   (*ti)=t_pre+t1;
   return(0);
}

/** 
  This function determines if the point A is exactly on the section.

  \param[in] sec 	type of Poincare section (sect = SEC1, SEC2, or SECg).
  \param[in] a 		point, 4 coordinates: (l, L, g, G).

  \returns true if point $a$ is exactly on section, false if it is not.
  */
bool onsection_del_car (section_t sec, double a[DIM])
{
   bool bonsection = false;
   double l,g;

   l = a[0];
   g = a[2];

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
      case SECg :       // section {g=0}
         {
            bonsection = (fabs(remainder(g,TWOPI))<POINCARE_DEL_CAR_TOL);
            break;
         }
      case SECg2 :       // section {g=\pi}
         {
            bonsection = (fabs(remainder(g-M_PI,TWOPI))<POINCARE_DEL_CAR_TOL);
            break;
         }
   }
   return(bonsection);
}

double inter_del_car_f(double t, void *p)
{
   int i, status;
   double pt[DIM];      // point in Cartesian coordinates

   struct inter_del_car_f_params *params = (struct inter_del_car_f_params *)p;
   double mu = (params->mu);
   section_t sec = (params->sec);
   pt[0] = (params->x); 
   pt[1] = (params->y); 
   pt[2] = (params->px);
   pt[3] = (params->py);

   // auxiliary variables
   double pt_del[DIM];  // point in Delaunay coordinates
   double l,g;
   double d;            // distance to section

   status=frtbp(mu,t,pt);	// flow(t,pt)
   if(status!=0)
   {
      fprintf(stderr, "inter_del_car_f: error computing flow");
      exit(EXIT_FAILURE);
   }

   cardel(pt,pt_del);
   l = pt_del[0];
   g = pt_del[2];

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
      case SECg :       // section {g=0}
         {
            d = remainder(g,TWOPI);
	        //fprintf(stderr,"DEBUG: t=%.15le, y=%.15le, px=%.15le, g=%.15le\n",
			//		t, pt[1], pt[2], g);
            break;
         }
      case SECg2 :       // section {g=\pi}
         {
            d = remainder(g-M_PI,TWOPI);
            break;
         }
   }
   return(d);
}

// name OF FUNCTION: inter_del_car
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
// "inter_del_car_f" for the intersection time "t".
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

int inter_del_car(double mu, section_t sec, double epsabs, 
        double x_del[DIM], double x_car[DIM], 
        double t0, double t1, double *t)
{
    int status;
    int iter = 0, max_iter = 100;
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    double f;
    gsl_function F;

	//double t_lo = t0;
	//double t_hi = t1;

    struct inter_del_car_f_params params = {mu, sec, x_car[0], x_car[1], x_car[2], x_car[3]};
    *t=0.0;
  
    F.function = &inter_del_car_f;
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
             "inter_del_car: solver iteration failed, gsl_errno=%d\n", \
             status);
           return(1);
        }
        *t = gsl_root_fsolver_root (s);

		// Unfortunately, cardel does not have good numerical accuracy, and
		// thus the residual |cardel(frtbp(t))| will never get very small even
		// though the root is computed with extreme precision |t-t^*|<10^{-15}.
		// Thus we stop as soon as forward error is small enough.
        f=inter_del_car_f(*t,&params);
        status = gsl_root_test_residual(f, epsabs);
		//t_lo = gsl_root_fsolver_x_lower(s);
		//t_hi = gsl_root_fsolver_x_upper(s);
        //status = gsl_root_test_interval(t_lo, t_hi, epsabs, 0);
      }
    while (status == GSL_CONTINUE && iter < max_iter);
    gsl_root_fsolver_free (s);

    if(iter>=max_iter && f>1.e-10)
    {
       fprintf(stderr, \
               "inter_del_car: maximum number of iterations reached\n");
       fprintf(stderr, \
               "inter_del_car: latest residual: %.15e\n",f);
       return(ERR_MAXITER_DEL_CAR);
    }

    // "*t" is the intersection time.
    if(frtbp(mu,*t,x_car))	// compute the intersection point "x"
    {
       fprintf(stderr, "inter_del_car: error computing intersection point\n");
       exit(EXIT_FAILURE);
    }
    cardel(x_car,x_del);
    return 0;
}
