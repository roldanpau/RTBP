// FILE:          trtbp.c
// TITLE:         Trajectory of Restricted Three Body Problem
// AUTHOR:        Pau Roldan
//                All code is my own except where credited to others.
// DATE:          September 21, 2009
//
// PURPOSE:
// This program computes a trajectory of the RTBP given an initial condition
// (t_0,x_0) and an integration time (positive or negative).
// It reads the mass parameter, initial condition, integration time, and step
// size from stdin. 
// The trajectory corresponding to the initial condition is integrated
// numerically. 
// Points in the trajectory are written to stdout, in the following format:
//    t_0 x_0
//    t_1 x_1
//    t_2 x_2
//    ...
//    t_N x_N,
//  where t_N=t_0+integration_time, and x_i is the position at time t_i.
//  Points in the trajectory are equispaced, meaning that t_{i+1} - t_i =
//  constant step size.
//
// NOTES:
// Here, step size does not refer to integration step size, which is actually
// variable, but to the distance between points in trajectory.
//
// The step size must have the same sign as the integration time!!!
// (positive->fwd integration, negative->bwd integration).
//
// OVERALL METHOD:
// The list of general tasks is:
// 1. Input mass parameter, initial condition, integration time, and step
//    size from stdin.
// 2. Integrate trajectory numerically by calling frtbp. 
//    As the trajectory is integrated, output points in the trajectory to
//    stdout.
//
// FUNCTIONS:
//
// print_pt
//    Print a point of a trajectory to stdout.

#include <stdio.h>
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE
#include <math.h>	// fabs
#include <frtbp.h>	// frtbp
#include <rtbp.h>	// DIM

// name OF FUNCTION: print_pt
// CREDIT: 
// PURPOSE:
// Print a point of a trajectory to stdout.
// The format in which the point is printed is:
//    t x\n
// where x is the position at time t.
//
// PARAMETERS:
// - t adimensional time at which the trajectory is at point x.
// - x point in phase space, 4 coordinates: (x, X, y, Y).
// 
// RETURN VALUE:
//
// CALLS TO: none
//
// CALLED FROM: main

void print_pt(double t, double x[DIM])
{
   if(printf("%le %le %le %le %le\n", t, x[0], x[1], x[2], x[3])<0)
   {
      perror("main: error writting output");
      exit(EXIT_FAILURE);
   }
}

int main( )
{
   double mu;		/* mass parameter */
   double t;		/* current time */
   double t0;		/* initial time */
   double t1;		/* final time */
   double h;		/* step size */
   double x[DIM];
   int status;

   // Input mass parameter, initial condition, integration time, and step
   // size from stdin.
   if(scanf("%le %le %le %le %le %le %le %le", &mu, &t0, x, x+1, x+2, x+3, 
	    &t1, &h) < 8)
   {
      perror("main: error reading input");
      exit(EXIT_FAILURE);
   }

   // The step size must have the same sign as the integration time!!!
   // (positive->fwd integration, negative->bwd integration).
   if(t1*h<0)
   {
      fprintf(stderr, \
	    "main: step size must have the same sign as integration time\n");
      exit(EXIT_FAILURE);
   }
	    
   // Integrate trajectory numerically.
   t=t0;
   do
   {
      // Output time and point to stdout.
      print_pt(t, x);

      status = frtbp(mu,h,x);
      if (status != 0)
      {
	 fprintf(stderr, "main: error integrating trajectory");
	 exit(EXIT_FAILURE);
      }
      t += h;
   } while (fabs(t-t0)<fabs(t1));

   // Final point is treated specially.
   status = frtbp(mu,-h,x);
   t -= h;
   status = frtbp(mu,t1-(t-t0),x);
   if (status != 0)
   {
      fprintf(stderr, "main: error integrating trajectory");
      exit(EXIT_FAILURE);
   }
   t = t0+t1;
   print_pt(t, x);
   exit(EXIT_SUCCESS);
}

