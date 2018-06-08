/*! \file
    \brief Restricted Three Body Problem equations in Delaunay coordinates

    $Author: roldan $
    $Date: 2013-03-26 22:26:03 $
*/

// FUNCTIONS
// =========
//
// Hamilt_del
//    Computes the Hamiltonian of the RTBP problem in Delaunay coordinates.
//
// re_DHell 
// im_DHell
//    This function computes the function $\Delta_{ell}^{1,+}$ (real and
//    imaginary part).

#include <math.h>		// M_PI
#include <stdlib.h>		// EXIT_FAILURE
#include <assert.h>
#include <gsl/gsl_errno.h>	// GSL_SUCCESS
#include <gsl/gsl_roots.h>	// GSL one-dimensional root finding
#include "rtbpdel.h"		// ERR_COLLISION

//const double COLLISION_TOL = 1.e-12;

struct trig_params
{
   double e, l;
};

double trig(double u, void *params)
{
   struct trig_params *p = (struct trig_params *)params;

   double e = p->e;
   double l = p->l;

   return l - u + e*sin(u);
}

/*
double trig_deriv(double u, void *params)
{
   struct trig_params *p = (struct trig_params *)params;

   double e = p->e;
   double l = p->l;

   return -1.0 + e*cos(u);
}

void trig_fdf(double u, void *params, double *y, double *dy)
{
   struct trig_params *p = (struct trig_params *)params;

   double e = p->e;
   double l = p->l;

   *y = l - u + e*sin(u);
   *dy = -1.0 + e*cos(u);
}
*/

// Compute the eccentric anomaly u.
// We need to solve u-e sin(u) = l.
// We don't impose that u is in [0,2pi) on output.  (If l<0 then u may be
// negative as well.)
double eccentric(double e, double l)
{
   // Desired precision for root. 
   // Note: we have tried 1.e-15, but such high precision could not always be
   // achieved.
   const double epsabs = 1.e-15;

    int status;
    int iter = 0, max_iter = 1000;
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    double r = 0;
    double x_lo = l-1, x_hi = 1+l;
    gsl_function F;
    struct trig_params params = {e, l};
  
    F.function = &trig;
    F.params = &params;
  
    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc (T);
    gsl_root_fsolver_set (s, &F, x_lo, x_hi);
  
    /*
    printf ("using %s method\n", 
	    gsl_root_fsolver_name (s));
  
    printf ("%5s [%9s, %9s] %9s %9s\n",
	    "iter", "lower", "upper", "root", 
	    "err(est)");
	    */
  
    do
      {
	iter++;
	status = gsl_root_fsolver_iterate (s);
	r = gsl_root_fsolver_root (s);
	x_lo = gsl_root_fsolver_x_lower (s);
	x_hi = gsl_root_fsolver_x_upper (s);
	status = gsl_root_test_interval (x_lo, x_hi,
					 epsabs, 0);
  
	/*
	if (status == GSL_SUCCESS)
	  printf ("Converged:\n");
  
	printf ("%5d [%.7f, %.7f] %.7f %.7f\n",
		iter, x_lo, x_hi,
		r,
		x_hi - x_lo);
		*/
      }
    while (status == GSL_CONTINUE && iter < max_iter);
  
    if(status != GSL_SUCCESS)
    {
       fprintf(stderr, "eccentric: can not find root!\n");
       exit(EXIT_FAILURE);
    }
    gsl_root_fsolver_free (s);

    //assert((0<=r) && (r<2*M_PI));
    return(r);
}

// This is used in functions re_DHell, im_DHell below
double Delta(double r, double v, double g)
{
   return sqrt(r*r + 1.0 - 2.0*r*cos(v+g));
}

double N_func(double r, double v, double g)
{
   return 1.0/sqrt(r*r + 1.0 - 2.0*r*cos(v+g));
}

double R(double mu, const double *x)
{
   double l = x[0];
   double L = x[1];
   double g = x[2];
   double G = x[3];

   // eccentricity
   double e = sqrt(1.0 - G*G/(L*L));

   // eccentric anomaly (angle u)
   double u = eccentric(e,l);

   // NOTE: v in (-pi,pi)
   double v = 2.0*atan( sqrt((1.0+e)/(1.0-e))*tan(u/2.0) );

   // modulus of asteroid r
   double r = L*L*(1.0-e*cos(u));

   // auxiliary variables
   double umu = 1.0-mu;

   return -umu/mu*N_func(-r/mu,v,g) - mu/umu*N_func(r/umu,v,g) + 1.0/r;
}

double Hamilt_del(double mu, const double *x)
{
   double L = x[1];
   double G = x[3];

   return -1.0/(2*L*L) - G + R(mu,x);
}

double dN_r(double r, double v, double g)
{
   double N = N_func(r,v,g);
   return (cos(v+g)-r)*N*N*N;
}

double dN_v(double r, double v, double g)
{
   double N = N_func(r,v,g);
   return -r*sin(v+g)*N*N*N;
}

int rtbp_del(double t, const double *x, double *y, void *params)
{
   double mu = *(double *)params;
   double l = fmod(x[0],2*M_PI);
   double L = x[1];
   double g = x[2];
   double G = x[3];

   // eccentricity
   double e = sqrt(1.0 - G*G/(L*L));

   // eccentric anomaly (angle u)
   double u = eccentric(e,l);

   // NOTE: v in (-pi,pi)
   double v = 2.0*atan( sqrt((1.0+e)/(1.0-e))*tan(u/2.0) );

   // auxiliary variables
   double su = sin(u);
   double cu = cos(u);
   double sv = sin(v);
   double cv = cos(v);

   // modulus of asteroid r
   double r = L*L*(1.0-e*cu);

   // auxiliary variables
   double umu = 1.0-mu;
   double musq = mu*mu;

   double c1 = umu/musq;
   double c2 = -mu/(umu*umu);

   double c3 = -umu/mu;
   double c4 = -mu/umu;
   double c5 = -1/(r*r);
   double Gsq = G*G;
   double esq = e*e;
   double Lsq = L*L;
   double Lcu = Lsq*L;

   // partial derivatives of N evaluated at -r/mu
   double dN_r_mu = dN_r(-r/mu,v,g);
   double dN_v_mu = dN_v(-r/mu,v,g);
   double dN_g_mu = dN_v_mu;

   // partial derivatives of N evaluated at r/(1-mu)
   double dN_r_umu = dN_r(r/umu,v,g);
   double dN_v_umu = dN_v(r/umu,v,g);
   double dN_g_umu = dN_v_umu;

   // partial derivatives of r
   double dr_G = G*cv/e;
   double dr_L = 1.0/L*(2.0*r-Gsq*cv/e);
   double dr_l = Lsq*e*sv/sqrt(1.0-esq);

   // partial derivatives of v
   double dv_G = -sv/(e*G)*(2.0+e*cv);
   double dv_L = sv/(1.0-esq)*(2.0+e*cv)*Gsq/(e*Lcu);
   double dv_l = sqrt(1.0-esq)/((1.0-e*cu)*(1.0-e*cu));

   // partial derivatives of R
   double dR_G = (c1*dN_r_mu + c2*dN_r_umu)*dr_G + 
      (c3*dN_v_mu + c4*dN_v_umu)*dv_G +
      c5*dr_G;
   double dR_L = (c1*dN_r_mu + c2*dN_r_umu)*dr_L + 
      (c3*dN_v_mu + c4*dN_v_umu)*dv_L +
      c5*dr_L;
   double dR_l = (c1*dN_r_mu + c2*dN_r_umu)*dr_l + 
      (c3*dN_v_mu + c4*dN_v_umu)*dv_l +
      c5*dr_l;
   double dR_g = c3*dN_g_mu + c4*dN_g_umu;

   // vector field
   y[0] = 1.0/Lcu + dR_L;	// \dot l
   y[1] = -dR_l;		// \dot L
   y[2] = -1.0 + dR_G;		// \dot g
   y[3] = -dR_g;		// \dot G

   // DEBUG:
   // fprintf(stderr, "%e %e %e\n", t, y[0], y[2]);
   
   return GSL_SUCCESS;
}

int dot_l(const double *x, double *dl, void *params)
{
   double mu = *(double *)params;
   double l = fmod(x[0],2*M_PI);
   double L = x[1];
   double g = x[2];
   double G = x[3];

   // eccentricity
   double e = sqrt(1.0 - G*G/(L*L));

   // eccentric anomaly (angle u)
   double u = eccentric(e,l);

   // NOTE: v in (-pi,pi)
   double v = 2.0*atan( sqrt((1.0+e)/(1.0-e))*tan(u/2.0) );

   // auxiliary variables
   double cu = cos(u);
   double sv = sin(v);
   double cv = cos(v);

   // modulus of asteroid r
   double r = L*L*(1.0-e*cu);

   // auxiliary variables
   double umu = 1.0-mu;
   double musq = mu*mu;

   double c1 = umu/musq;
   double c2 = -mu/(umu*umu);

   double c3 = -umu/mu;
   double c4 = -mu/umu;
   double c5 = -1/(r*r);
   double Gsq = G*G;
   double esq = e*e;
   double Lsq = L*L;
   double Lcu = Lsq*L;

   // partial derivatives of N evaluated at -r/mu
   double dN_r_mu = dN_r(-r/mu,v,g);
   double dN_v_mu = dN_v(-r/mu,v,g);

   // partial derivatives of N evaluated at r/(1-mu)
   double dN_r_umu = dN_r(r/umu,v,g);
   double dN_v_umu = dN_v(r/umu,v,g);

   // partial derivatives of r
   double dr_L = 1.0/L*(2.0*r-Gsq*cv/e);

   // partial derivatives of v
   double dv_L = sv/(1.0-esq)*(2.0+e*cv)*Gsq/(e*Lcu);

   // partial derivatives of R
   double dR_L = (c1*dN_r_mu + c2*dN_r_umu)*dr_L + 
      (c3*dN_v_mu + c4*dN_v_umu)*dv_L +
      c5*dr_L;

   // vector field
   *dl = 1.0/Lcu + dR_L;	// \dot l

   return GSL_SUCCESS;
}

int dot_g(const double *x, double *dg, void *params)
{
   double mu = *(double *)params;
   double l = fmod(x[0],2*M_PI);
   double L = x[1];
   double g = x[2];
   double G = x[3];

   // eccentricity
   double e = sqrt(1.0 - G*G/(L*L));

   // eccentric anomaly (angle u)
   double u = eccentric(e,l);

   // NOTE: v in (-pi,pi)
   double v = 2.0*atan( sqrt((1.0+e)/(1.0-e))*tan(u/2.0) );

   // auxiliary variables
   double cu = cos(u);
   double sv = sin(v);
   double cv = cos(v);

   // modulus of asteroid r
   double r = L*L*(1.0-e*cu);

   // auxiliary variables
   double umu = 1.0-mu;
   double musq = mu*mu;
   double c1 = umu/musq;
   double c2 = -mu/(umu*umu);
   double c3 = -umu/mu;
   double c4 = -mu/umu;
   double c5 = -1/(r*r);

   // partial derivatives of N evaluated at -r/mu
   double dN_r_mu = dN_r(-r/mu,v,g);
   double dN_v_mu = dN_v(-r/mu,v,g);

   // partial derivatives of N evaluated at r/(1-mu)
   double dN_r_umu = dN_r(r/umu,v,g);
   double dN_v_umu = dN_v(r/umu,v,g);

   // partial derivatives of r
   double dr_G = G*cv/e;

   // partial derivatives of v
   double dv_G = -sv/(e*G)*(2.0+e*cv);

   // partial derivatives of R
   double dR_G = (c1*dN_r_mu + c2*dN_r_umu)*dr_G +
      (c3*dN_v_mu + c4*dN_v_umu)*dv_G +
      c5*dr_G;

   // vector field
   *dg = -1.0 + dR_G;		// \dot g

   return GSL_SUCCESS;
}

// Compute \mu\partial_G \Delta H_{\circ}.
// This function is almost identical to dot_g, except that we avoid the
// addition -1 + \mu\partial_G \Delta H_{\circ}, since we would be loosing
// accuracy (1 and \mu\partial_G \Delta H_{\circ} differ by orders of
// magnitude).
int mu_dDHcirc_G(const double *x, double *res, void *params)
{
   double mu = *(double *)params;
   double l = fmod(x[0],2*M_PI);
   double L = x[1];
   double g = x[2];
   double G = x[3];

   // eccentricity
   double e = sqrt(1.0 - G*G/(L*L));

   // eccentric anomaly (angle u)
   double u = eccentric(e,l);

   // NOTE: v in (-pi,pi)
   double v = 2.0*atan( sqrt((1.0+e)/(1.0-e))*tan(u/2.0) );

   // auxiliary variables
   double cu = cos(u);
   double sv = sin(v);
   double cv = cos(v);

   // modulus of asteroid r
   double r = L*L*(1.0-e*cu);

   // auxiliary variables
   double umu = 1.0-mu;
   double musq = mu*mu;
   double c1 = umu/musq;
   double c2 = -mu/(umu*umu);
   double c3 = -umu/mu;
   double c4 = -mu/umu;
   double c5 = -1/(r*r);

   // partial derivatives of N evaluated at -r/mu
   double dN_r_mu = dN_r(-r/mu,v,g);
   double dN_v_mu = dN_v(-r/mu,v,g);

   // partial derivatives of N evaluated at r/(1-mu)
   double dN_r_umu = dN_r(r/umu,v,g);
   double dN_v_umu = dN_v(r/umu,v,g);

   // partial derivatives of r
   double dr_G = G*cv/e;

   // partial derivatives of v
   double dv_G = -sv/(e*G)*(2.0+e*cv);

   // partial derivatives of R
   double dR_G = (c1*dN_r_mu + c2*dN_r_umu)*dr_G +
      (c3*dN_v_mu + c4*dN_v_umu)*dv_G +
      c5*dr_G;

   // Result
   *res = dR_G;		// \mu \partial_G \Delta H_{\circ}

   return GSL_SUCCESS;
}

// NOTES:
//    We follow the convention to place the large mass (Sun) to the left of
//    the origin, and the small mass (Jupiter) to the right.
//    This is opposite to the usual astrodynamics convention.
//
//    This function is used by program inner_circ, but we put it here because
//    it is very similar to function rtbp_del.
//
//    We normalize $\ell$ between 0 and 2\pi to compute the vectorfield.
//    This is done so that function "eccentric" is more precise.
//
// CALLED FROM: outer_circ::integrand_omega_in

int f0(const double *x, double *res, void *params)
{
   // auxiliary variables
   double dl;	// \dot l

   dot_l(x,&dl,params);

   // recall that f0 = (3-l'(x))/(3 l'(x))
   *res = (3.0-dl)/(3.0*dl);		// f0

   return GSL_SUCCESS;
}

// NOTES:
//    We follow the convention to place the large mass (Sun) to the left of
//    the origin, and the small mass (Jupiter) to the right.
//    This is opposite to the usual astrodynamics convention.
//
//    This function is used by program inner_circ, but we put it here because
//    it is very similar to function rtbp_del.
//
// Careful!!! 
// Marcel's formulas factorize \mu in front of the integral.  We include \mu
// in the integrand.
//
//    We normalize $\ell$ between 0 and 2\pi to compute the vectorfield.
//    This is done so that function "eccentric" is more precise.
//
// CALLED FROM: outer_circ_stoch::integrand_omega_pm

int f0_stoch(const double *x, double *res, void *params)
{
   // auxiliary variables
   double mu_dDH_G;	// \mu \partial_G \Delta H_{\circ}

   mu_dDHcirc_G(x,&mu_dDH_G,params);

   // recall that f0 = \mu\partial_G \Delta H_{\circ}(x) / 
   //   (-1+\mu\partial_G \Delta H_{\circ})(x).
   *res = mu_dDH_G/(-1+mu_dDH_G);		// f0

   return GSL_SUCCESS;
}

// name OF FUNCTION: re_DHell
// CREDIT: Marcel Guardia and Pau Roldan
// PURPOSE:
//    The inner dynamics of the elliptic problem is given by the complex
//    integral $A^+$ (see Marcel's notes "Inner and outer dynamics").
//    This function computes the function $\Delta H_{ell}^{1,+}$ (real part).
//
// NOTES:
//    This function is located in this file because it closely resembles the
//    function rtbp_del.
//
//    We normalize $\ell$ between 0 and 2\pi to compute the vectorfield.
//    This is done so that function "eccentric" is more precise.
//
// PARAMETERS:
// - x point in phase space, 4 coordinates: (l,L,g,G).
// - params pointer to the parameter of the system: the mass ratio "mu".
// 
// RETURN VALUE:
// value of the function $\Delta H_{ell}^{1,+}$ (real part).
//
// CALLS TO:
//
// CALLED FROM: re_integrand_inner_ell

double re_DHell(const double *x, void *params)
{
   double mu = *(double *)params;
   double l = fmod(x[0],2*M_PI);
   double L = x[1];
   double g = x[2];
   double G = x[3];

   // eccentricity
   double e = sqrt(1.0 - G*G/(L*L));

   // eccentric anomaly (angle u)
   double u = eccentric(e,l);

   // NOTE: v in (-pi,pi)
   double v = 2.0*atan( sqrt((1.0+e)/(1.0-e))*tan(u/2.0) );

   // modulus of asteroid r
   double r = L*L*(1.0-e*cos(u));

   // auxiliary variables
   double umu = 1.0-mu;

   // auxiliary variables
   double vg = v+g;

   double D_r_mu = Delta(-r/mu,v,g);
   double D_r_umu = Delta(r/umu,v,g);

   double D_r_mu_cb = D_r_mu*D_r_mu*D_r_mu;
   double D_r_umu_cb = D_r_umu*D_r_umu*D_r_umu;

   return umu/mu*(1+r/mu*cos(vg))/(2*D_r_mu_cb) + \
       mu/umu*(1-r/umu*cos(vg))/(2*D_r_umu_cb);
}

// name OF FUNCTION: im_DHell
// CREDIT: Marcel Guardia and Pau Roldan
// PURPOSE:
//    The inner dynamics of the elliptic problem is given by the complex
//    integral $A^+$ (see Marcel's notes "Inner and outer dynamics").
//    This function computes the function $\Delta H_{ell}^{1,+}$
//    (imarinary part). See also my yellow notebook.
//
// NOTES:
//    This function is located in this file because it closely resembles the
//    function rtbp_del.
//
//    We normalize $\ell$ between 0 and 2\pi to compute the vectorfield.
//    This is done so that function "eccentric" is more precise.
//
// PARAMETERS:
// - x point in phase space, 4 coordinates: (l,L,g,G).
// - params pointer to the parameter of the system: the mass ratio "mu".
// 
// RETURN VALUE:
// value of the function $\Delta H_{ell}^{1,+}$ (imaginary part).
//
// CALLS TO:
//
// CALLED FROM: im_integrand_inner_ell

double im_DHell(const double *x, void *params)
{
   double mu = *(double *)params;
   double l = fmod(x[0],2*M_PI);
   double L = x[1];
   double g = x[2];
   double G = x[3];

   // eccentricity
   double e = sqrt(1.0 - G*G/(L*L));

   // eccentric anomaly (angle u)
   double u = eccentric(e,l);

   // NOTE: v in (-pi,pi)
   double v = 2.0*atan( sqrt((1.0+e)/(1.0-e))*tan(u/2.0) );

   // modulus of asteroid r
   double r = L*L*(1.0-e*cos(u));

   // auxiliary variables
   double umu = 1.0-mu;

   // auxiliary variables
   double vg = v+g;

   double D_r_mu = Delta(-r/mu,v,g);
   double D_r_umu = Delta(r/umu,v,g);

   double D_r_mu_cb = D_r_mu*D_r_mu*D_r_mu;
   double D_r_umu_cb = D_r_umu*D_r_umu*D_r_umu;

   return umu/mu*r/mu*sin(vg)/D_r_mu_cb - mu/umu*r/umu*sin(vg)/D_r_umu_cb;
}

// name OF FUNCTION: re_dDHell
// CREDIT: Marcel Guardia and Pau Roldan
// PURPOSE:
//    The inner dynamics of the elliptic problem is given by the complex
//    integral $A^+$ (see Marcel's notes "Inner and outer dynamics").
//    This function computes the function $\partial_t \Delta H_{ell}^{1,+}$
//    (real part). See also my yellow notebook.
//
// NOTES:
//    This function is located in this file because it closely resembles the
//    function rtbp_del.
//
//    We normalize $\ell$ between 0 and 2\pi to compute the vectorfield.
//    This is done so that function "eccentric" is more precise.
//
// PARAMETERS:
// - x point in phase space, 4 coordinates: (l,L,g,G).
// - params pointer to the parameter of the system: the mass ratio "mu".
// 
// RETURN VALUE:
// value of the function $\partial_t \Delta H_{ell}^{1,+}$ (real part).
//
// CALLS TO:
//
// CALLED FROM: re_integrand_inner_ell

/* This function should have never been used. Use re_DHell instead!
double re_dDHell(const double *x, void *params)
{
   double mu = *(double *)params;
   double l = fmod(x[0],2*M_PI);
   double L = x[1];
   double g = x[2];
   double G = x[3];

   // eccentricity
   double e = sqrt(1.0 - G*G/(L*L));

   // eccentric anomaly (angle u)
   double u = eccentric(e,l);

   // NOTE: v in (-pi,pi)
   double v = 2.0*atan( sqrt((1.0+e)/(1.0-e))*tan(u/2.0) );

   // modulus of asteroid r
   double r = L*L*(1.0-e*cos(u));

   // auxiliary variables
   double umu = 1.0-mu;
   double musq = mu*mu;

   // auxiliary variables
   double vg = v+g;

   double D_r_mu = Delta(-r/mu,v,g);
   double D_r_umu = Delta(r/umu,v,g);

   return -umu/musq*(-r/mu*sin(vg))/(D_r_mu*D_r_mu*D_r_mu) -
      mu/umu*(r/umu*sin(vg))/(D_r_umu*D_r_umu*D_r_umu);
}
*/

// name OF FUNCTION: im_dDHell
// CREDIT: Marcel Guardia and Pau Roldan
// PURPOSE:
//    The inner dynamics of the elliptic problem is given by the complex
//    integral $A^+$ (see Marcel's notes "Inner and outer dynamics").
//    This function computes the function $\partial_t \Delta H_{ell}^{1,+}$
//    (imarinary part). See also my yellow notebook.
//
// NOTES:
//    This function is located in this file because it closely resembles the
//    function rtbp_del.
//
//    We normalize $\ell$ between 0 and 2\pi to compute the vectorfield.
//    This is done so that function "eccentric" is more precise.
//
// PARAMETERS:
// - x point in phase space, 4 coordinates: (l,L,g,G).
// - params pointer to the parameter of the system: the mass ratio "mu".
// 
// RETURN VALUE:
// value of the function $\partial_t \Delta H_{ell}^{1,+}$ (imaginary part).
//
// CALLS TO:
//
// CALLED FROM: im_integrand_inner_ell

/* This function should have never been used. Use im_DHell instead!
double im_dDHell(const double *x, void *params)
{
   double mu = *(double *)params;
   double l = fmod(x[0],2*M_PI);
   double L = x[1];
   double g = x[2];
   double G = x[3];

   // eccentricity
   double e = sqrt(1.0 - G*G/(L*L));

   // eccentric anomaly (angle u)
   double u = eccentric(e,l);

   // NOTE: v in (-pi,pi)
   double v = 2.0*atan( sqrt((1.0+e)/(1.0-e))*tan(u/2.0) );

   // modulus of asteroid r
   double r = L*L*(1.0-e*cos(u));

   // auxiliary variables
   double umu = 1.0-mu;
   double musq = mu*mu;

   // auxiliary variables
   double vg = v+g;

   double D_r_mu = Delta(-r/mu,v,g);
   double D_r_umu = Delta(r/umu,v,g);

   return -umu/musq*(1.0+r/mu*cos(vg))/(2.0*D_r_mu*D_r_mu*D_r_mu) -
      mu/umu*(1.0-r/umu*cos(vg))/(2.0*D_r_umu*D_r_umu*D_r_umu);
}
*/
