// =======================================
// Restricted Three Body Problem equations
// =======================================
// FILE:          $RCSfile: rtbp.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2012-11-26 10:48:15 $
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
// Hamilt
//    Hamiltonian function of the RTBP.
//
// rtbp
//    Computes the vectorfield of the RTBP problem.
//
// rtbp_inv
//    Computes the negative vectorfield of the RTBP problem.

#include <math.h>
#include <gsl/gsl_errno.h>
#include "rtbp.h"		// ERR_COLLISION

const double COLLISION_TOL = 1.e-12;	

double Hamilt(double mu, const double *p)
{
   double x=p[0];
   double y=p[1];
   double px=p[2];
   double py=p[3];

   // Place large mass to the left of the origin, small mass to the right.
   double mu1 = mu;
   double mu2 = 1.0-mu;

   double r1=sqrt((x-mu2)*(x-mu2)+y*y);
   double r2=sqrt((x+mu1)*(x+mu1)+y*y);

   return 0.5*(px*px + py*py) + y*px - x*py - mu1/r1 - mu2/r2;
}

// name OF FUNCTION: rtbp
// CREDIT: Angel Jorba, with modifications by Pau Roldan
// PURPOSE:
//    The function computes the vectorfield of the RTBP at a point named "x".
//    Specifically, we use the equations of motion of the spatial circular
//    RTBP, in the rotating (synodic) coordinate system, using adimensional
//    quantities, and Hamiltonian formulation (position-momentum).
//    See Szebehelly, page 349. 
//    See also the notes on normal forms by Maciej Capinski.
//    The function sets the value of the vectorfield in the variable "y".
//    It returns a status code (success/error).
//
// NOTES
// -----
//
// We follow the convention to place the large mass (Sun) to the LEFT of the
// origin, and the small mass (Jupiter) to the RIGHT.
//
// Collision tolerance: we consider a collision if 
// r1^3 < COLLISION_TOL or r2^3 < COLLISION_TOL.
//
// PARAMETERS:
// - t adimensional time at which the vectorfield is evaluated. Since this is
//   an autonomous ODE (does not depend on time), this parameter is not used.
// - x point in phase space, 4 coordinates: (X, Y, P_X, P_Y).
// - y vectorfield at (t,x), 4 coordinates: d/dt(X, Y, P_X, P_Y).
// - params pointer to the parameter of the system: the mass ratio "mu".
// 
// RETURN VALUE:
// status code of the function (success/error):
//    - GSL_SUCCESS: success.
//    - ERR_COLLISION: collision of the third mass with one of the primaries.
//
// CALLS TO: none
//
// CALLED FROM: main

int rtbp(double t, const double *x, double *y, void *params)
{
   double r1,r13,r2,r23,aux3;
   double mu = *(double *)params;
   double mu1 = mu;
   double mu2 = 1.0-mu;
   r1=sqrt((x[0]-mu2)*(x[0]-mu2)+x[1]*x[1]);
   r13=r1*r1*r1;
   r2=sqrt((x[0]+mu1)*(x[0]+mu1)+x[1]*x[1]);
   r23=r2*r2*r2;
   if(r13<COLLISION_TOL || r23<COLLISION_TOL)
      return ERR_COLLISION;
   aux3=mu1/r13+mu2/r23;
   y[0]=x[2]+x[1];
   y[1]=-x[0]+x[3];
   y[2]=x[3]-mu1*(x[0]-mu2)/r13-mu2*(x[0]+mu1)/r23;
   y[3]=-x[2]-aux3*x[1];
   return GSL_SUCCESS;
}

// name OF FUNCTION: rtbp_inv
// CREDIT: Angel Jorba, with modifications by Pau Roldan
// PURPOSE:
//    The function computes the vectorfield of the RTBP at a point named "x".
//    Specifically, we use the equations of motion of the spatial circular
//    RTBP, in the rotating (synodic) coordinate system, using adimensional
//    quantities, and Hamiltonian formulation (position-momentum).
//    See Szebehelly, page 349. 
//    See also the notes on normal forms by Maciej Capinski.
//    The function sets the value of the vectorfield in the variable "y".
//    It returns a status code (success/error).
//
// NOTES
// -----
//
// We follow the convention to place the large mass (Sun) to the LEFT of the
// origin, and the small mass (Jupiter) to the RIGHT.
//
// Collision tolerance: we consider a collision if 
// r1^3 < COLLISION_TOL or r2^3 < COLLISION_TOL.
//
// PARAMETERS:
// - t adimensional time at which the vectorfield is evaluated. Since this is
//   an autonomous ODE (does not depend on time), this parameter is not used.
// - x point in phase space, 4 coordinates: (X, Y, P_X, P_Y).
// - y vectorfield at (t,x), 4 coordinates: d/dt(X, Y, P_X, P_Y).
// - params pointer to the parameter of the system: the mass ratio "mu".
// 
// RETURN VALUE:
// status code of the function (success/error):
//    - GSL_SUCCESS: success.
//    - ERR_COLLISION: collision of the third mass with one of the primaries.
//
// CALLS TO: none
//
// CALLED FROM: main

int rtbp_inv(double t, const double *x, double *y, void *params)
{
   double r1,r13,r2,r23,aux3;
   double mu = *(double *)params;
   double mu1 = mu;
   double mu2 = 1.0-mu;
   r1=sqrt((x[0]-mu2)*(x[0]-mu2)+x[1]*x[1]);
   r13=r1*r1*r1;
   r2=sqrt((x[0]+mu1)*(x[0]+mu1)+x[1]*x[1]);
   r23=r2*r2*r2;
   if(r13<COLLISION_TOL || r23<COLLISION_TOL)
      return ERR_COLLISION;
   aux3=mu1/r13+mu2/r23;
   y[0]=-x[2]-x[1];
   y[1]=x[0]-x[3];
   y[2]=-x[3]+mu1*(x[0]-mu2)/r13+mu2*(x[0]+mu1)/r23;
   y[3]=x[2]+aux3*x[1];
   return GSL_SUCCESS;
}
