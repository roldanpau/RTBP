// ===================================================
// Poincare map 2 of the Restricted Three Body Problem
// ===================================================
// FILE:          $RCSfile: prtbp2.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2012-06-06 13:03:01 $
//
// PURPOSE:
//
// NOTES:
//
// OVERALL METHOD:
//
// FUNCTIONS
// =========
//
// prtbp2
//    Compute the n-th iterate of the Poincare map $P^n(x)$ of the RTBP.
// prtbp2_inv
//    Compute the n-th iterate $P^{-n}(x)$ of the inverse Poincare map of the
//    RTBP.


#include <stdio.h>	// fprintf
#include <stdlib.h>	// EXIT_FAILURE
#include <gsl/gsl_errno.h>	// GSL_SUCCESS
#include <gsl/gsl_roots.h>
#include <frtbp.h>	// frtbp
#include <rtbp.h>	// DIM
#include <prtbp.h>	// POINCARE_TOL, inter

// name OF FUNCTION: prtbp2
// CREDIT: 
//
// PURPOSE
// =======
// Consider the RTBP in rotating coordinates.
// Let S be the Poincare section {y=0, p_y<0}.
// Let $x$ be a point in the section, and suppose that flow at a point $x$ is
// transversal to the section.
// Compute the n-th iterate of the Poincare map $P^n(x)$ of the RTBP.
// This procedure also computes, as a side product, the integration time to
// intersect the Poincare section "n" times.
// 
// PARAMETERS
// ==========
// mu
//    mass parameter for the RTBP
// cuts
//    number of iterates of the Poincare map (cuts=n: n cuts with the
//    Poincare section).
// x
//    Initial point, 4 coordinates: (X, Y, P_X, P_Y). 
//    On return of the this function, it holds the image point $P(x)$.
// ti
//    On return, it holds the integration time to intersect the Poincare
//    section "n" times.
// 
// RETURN VALUE
// ============
// Returns a non-zero error code to indicate an error and 0 to indicate
// success.
// If an integration error is encountered, the function returns a non-zero
// value.
//
// NOTES
// =====
// We do not impose that $x$ is on the section.
// We do not check that flow at $x$ is transversal to section!!!
//
// A point is assumed to be on the Poincare section if it is within distance
// POINCARE_TOL to the section.
//
// On successful return of this function, we set coordinate $y$ exactly equal
// to zero.
//
// CALLS TO: frtbp, inter

int prtbp2(double mu, int cuts, double x[DIM], double *ti)
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
      // Integrate trajectory until it reaches section y=0
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
	 status = frtbp(mu,0.01,x);
	 t += 0.01;
	 if (status != GSL_SUCCESS)
	 {
	    fprintf(stderr, "prtbp2: error integrating trajectory\n");
	    return(1);
	 }
      } 
      // while(no crossing of Poincare section), i.e.
      // while(!(y==0 && p_y<0))
      while(!((x[1]==0 || x_pre[1]*x[1]<0) && x[3]<0)); 
      n++;
   }
   if(x[1]==0)	// point "x" is exactly on the section {y=0}
   {
      (*ti)=t;
      return(0);
   }
   // Crossing happened between times t_pre and t. 

   // Restore previous value of point "x"
   for(i=0;i<DIM;i++)
      x[i]=x_pre[i];

   // Intersect trajectory starting at point x with section.
   if(inter(mu, POINCARE_TOL, x, 0, t-t_pre, &t1))
   {
      fprintf(stderr, "prtbp2: error intersectig trajectory with section\n");
      return(1);
   }
   // Here, coordinate y is 0 with tolerance POINCARE_TOL_DEL. 
   // We force y to be exactly equal to zero.
   x[1] = 0;    // y

   // Set time to reach Poincare section
   (*ti)=t_pre+t1;
   return(0);
}

// name OF FUNCTION: prtbp2_inv
// CREDIT: 
//
// PURPOSE
// =======
// Consider the RTBP in rotating coordinates.
// Let S be the Poincare section {y=0, p_y<0}.
// Let $x$ be a point in the section, and suppose that flow at a point $x$ is
// transversal to the section.
// Compute the n-th iterate $P^{-n}(x)$ of the inverse Poincare map of the
// RTBP.
// This procedure also computes, as a side product, the integration time to
// intersect the Poincare section "n" times.
// 
// PARAMETERS
// ==========
// mu
//    mass parameter for the RTBP
// cuts
//    number of iterates of the Poincare map (cuts=n: n cuts with the
//    Poincare section).
// x
//    Initial point, 4 coordinates: (X, Y, P_X, P_Y). 
//    On return of the this function, it holds the image point $P^{-1}(x)$.
// ti
//    On return, it holds the integration time to intersect the Poincare
//    section. Since we are computing the inverse Poincare map, this
//    integration time "ti" must be negative.
// 
// RETURN VALUE
// ============
// Returns a non-zero error code to indicate an error and 0 to indicate
// success.
// If an integration error is encountered, the function returns a non-zero
// value.
//
// NOTES
// =====
// We do not impose that $x$ is on the section.
// We do not check that flow at $x$ is transversal to section!!!
//
// A point is assumed to be on the Poincare section if it is within distance
// POINCARE_TOL to the section.
//
// On successful return of this function, we set coordinate $y$ exactly equal
// to zero.
//
// CALLS TO: frtbp, inter

int prtbp2_inv(double mu, int cuts, double x[DIM], double *ti)
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
      // Integrate trajectory backwards until it reaches section y=0
      do
      {
	 // Save previous value of point "x" and time "t"
	 for(i=0;i<DIM;i++)
	    x_pre[i]=x[i];
	 t_pre=t;

	 // Integrate backwards for a "short" time t1=-1.0, short enough so
	 // that we can detect crossing of Poincare section.
	 status = frtbp(mu,-0.01,x);
	 t -= 0.01;
	 if (status != GSL_SUCCESS)
	 {
	    fprintf(stderr, "prtbp2_inv: error integrating trajectory\n");
	    return(1);
	 }
      } 
      // while(no crossing of Poincare section), i.e.
      // while(!(y==0 && p_y<0))
      while(!((x[1]==0 || x_pre[1]*x[1]<0) && x[3]<0)); 
      n++;
   }
   if(x[1]==0)	// point "x" is exactly on the section {y=0}
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
   if(inter(mu, POINCARE_TOL, x, t-t_pre, 0, &t1))
   {
      fprintf(stderr, "prtbp2_inv: error intersectig trajectory with section\n");
      return(1);
   }
   // Here, coordinate y is 0 with tolerance POINCARE_TOL_DEL. 
   // We force y to be exactly equal to zero.
   x[1] = 0;    // y

   // Set time to reach Poincare section
   (*ti)=t_pre+t1;
   return(0);
}
