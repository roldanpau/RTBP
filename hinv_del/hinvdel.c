// ===========================
// Invert Hamiltonian equation
// ===========================
// FILE:          $RCSfile: hinvdel.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2012-06-06 10:56:19 $
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
// hinv
//    Invert Hamiltonian equation
//       H(l,L,g,G)=H0,
//    solving for the unknown L.

#include <stdio.h>	// fprintf
#include <math.h>	// pow
#include <rtbp.h>	// DIM
#include <rtbpdel.h>	// Hamilt_del, rtbp_del
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>	// root-finding routines

const double EPSABS_NEWTON = 1e-14;
const int MAXITER_NEWTON = 100;

struct DeltaH_params
{
   double mu, H0, l, g, G;
};

// name OF FUNCTION: DeltaH, DeltaH_deriv, DeltaH_fdf
//
// PURPOSE
// =======
// Given a point $(l,L,g,G)$ and an energy value $H_0$, compute the
// difference function
// \[ H(l,L,g,G) - H_0 \]
// and the derivative of this function with respect to $L$.
//
// NOTES
// =====
// This function is used in the root-finding algorithm in hinv_del below.
//
// PARAMETERS
// ==========
// L
//    L coordinate
// params
//    other parameters: mu, H0, l,g,G
// 
// RETURN VALUE
// ============
// Returns the difference $H(l,L,g,G) - H_0$.
//
// CALLS TO: Hamilt_del, rtbp_del
//
// CALLED FROM: hinv_del


double DeltaH(double L, void *params)
{
   struct DeltaH_params *p = (struct DeltaH_params *)params;
   double mu = p->mu;
   double H0 = p->H0;

   double l = p->l;
   double g = p->g;
   double G = p->G;

   double x[DIM] = {l,L,g,G};
   return Hamilt_del(mu,x) - H0;
}

double DeltaH_deriv(double L, void *params)
{
   struct DeltaH_params *p = (struct DeltaH_params *)params;
   double mu = p->mu;
   double H0 = p->H0;

   double l = p->l;
   double g = p->g;
   double G = p->G;

   double x[DIM] = {l,L,g,G};

   double y[DIM]; 	// vector field evaluated at x
   int status;

   // Note: for the moment, we do not code a separate function dot_l to
   // compute $\partial H / \partial L$. Instead, we call rtbp_del to compute
   // the whole vector field, and return the first component.
   status = rtbp_del(0,x,y,&mu);
   if(status != GSL_SUCCESS)
   {
      fprintf(stderr, "DeltaH_deriv: error computing vector field");
      exit(EXIT_FAILURE);
   }
   return y[0];
}

void DeltaH_fdf(double L, void *params, double *y, double *dy)
{
   *y = DeltaH(L,params);
   *dy = DeltaH_deriv(L,params);
}

// name OF FUNCTION: hinv_del
// CREDIT: 
//
// PURPOSE
// =======
// Consider the Hamiltonian $H$ of the RTBP in rotating Delaunay coordinates
// \[ H(l,L,g,G) = -1/(2L^2) -G + \mu \Delta H_{circ}. \]
// Let the value of the Hamiltonian H=H0 be given. Supose we know the value
// of three coordinates, l,g,G.
// This procedure inverts the Hamiltonian equation
//    H(l,L,g,G)=H0,
// solving for the unknown L.
//
// NOTES
// =====
// Since the Hamiltonian equation is not solvable analytically for L (unlike
// Cartesian coordinates), we solve it using a Newton method. 
// For the 3:1 resonance, L is close to 3^{-1/3} so we use this value as a
// first guess for the Newton method.
//
// PARAMETERS
// ==========
// mu
//    mass parameter for the RTBP
// H
//    energy value
// p
//    point, 4 coordinates: p=(l,L,g,G). 
//    On input, p[0],p[2],p[3] holds the known values of l, g and G.
//    On return of the this function, p[1] holds the value of L such that
//    H(l,L,g,G)=H.
// 
// RETURN VALUE
// ============
// Returns GSL_SUCCESS) to indicate success, and a non-zero error code
// otherwise.
// If an integration error is encountered, the function returns a non-zero
// value, and p is unmodified.

int hinv_del(double mu, double H,double p[DIM])
{
    int status;
    int iter = 0;
    const gsl_root_fdfsolver_type *T;
    gsl_root_fdfsolver *s;
    double L = pow(3.0, -1.0/3.0);
    gsl_function_fdf FDF;
    struct DeltaH_params params = {mu,H,p[0],p[2],p[3]};	// mu,H,l,g,G

    // auxiliary variables
    double residual;
  
    FDF.f = &DeltaH;
    FDF.df = &DeltaH_deriv;
    FDF.fdf = &DeltaH_fdf;
    FDF.params = &params;
  
    T = gsl_root_fdfsolver_newton;
    s = gsl_root_fdfsolver_alloc (T);
    gsl_root_fdfsolver_set (s, &FDF, L);
  
    //printf ("%-5s %10s %10s\n", "iter", "root", "err(est)");
    do
      {
	iter++;
	status = gsl_root_fdfsolver_iterate (s);
	L = gsl_root_fdfsolver_root (s);
	residual = DeltaH(L,&params);
	status = gsl_root_test_residual(residual, EPSABS_NEWTON);
  
	//if (status == GSL_SUCCESS)
	//  printf ("Converged:\n");
  
	//printf ("%5d %10.7f %10.7f\n", iter, L, residual);
      }
    while (status == GSL_CONTINUE && iter < MAXITER_NEWTON);
  
    if(status == GSL_SUCCESS) 
       p[1] = L;

    gsl_root_fdfsolver_free (s);
    return status;
}
