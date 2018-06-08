// ====================================================
// Outer Map of the Elliptic Problem (Stochastic paper)
// ====================================================
// FILE:          $RCSfile: outer_ell_stoch.c,v $
//
// PURPOSE
// -------
// This program computes the shift of the outer map of the elliptic problem,
// given by the integral
// 
// \[ B_out^j(I) = -\mu B^+(I) +\mu C^+(I), \]
//
// where 
// 
// \[ B^+(I) = i\int_0^{2M\pi} 
//    f(\gamma_h(s)) e^{it(s)} - f(\gamma_p(s)) e^{i(t(s)+\omega_+^j)} ds, \]
// \[ C^+(I) = i\int_0^{-2M\pi}
//    f(\gamma_h(s)) e^{it(s)} - f(\gamma_p(s)) e^{i(t(s)+\omega_-^j)} ds. \]
//
// Here, $\gamma_h(s)$ is the homoclinic trajectory, and $\gamma_p(s)$ is the
// periodic trajectory.
// (see notes "Inner and outer dynamics" by Marcel).
//
// FUNCTIONS
// =========
//
// re_integrand_B_stoch
// im_integrand_B_stoch
// --------------
// Given the periodic point $p$ and the homoclinic point $h$ of the circular
// RTBP, an integration time $s$ in the homoclinic and periodic trajectory, 
// this function computes the integrand of the function $B^+$ or $C^+$
// (depending if this function is called from re_B_stoch or re_C_stoch),
// \[ f(\gamma_h(s)) e^{it(s)} - f(\gamma_p(s)) e^{i(t(s)+\omega_+^j)} \]
// \[ f(\gamma_h(s)) e^{it(s)} - f(\gamma_p(s)) e^{i(t(s)+\omega_-^j)} \]
//
// re_B_stoch
// im_B_stoch
// ----
// Given an energy level $H$, compute the integral $B^+$.
// This is computed using numerical integration.
//
// re_C_stoch
// im_C_stoch
// ----
// Given an energy level $H$, compute the integral $C^+$.
// This is computed using numerical integration.

#include <stdio.h>	// perror
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE

#include <gsl/gsl_integration.h>	// gsl_integration_qags

#include <utils_module.h>       // dblcpy
#include <rtbpdel.h>			// rtbp_del
#include <frtbpred.h>
#include <inner_ell_stoch.h>	// re_f_integrand_stoch, im_f_integrand_stoch

// 1.e-6 is too much
const double RELERROR = 1.e-2;
const size_t NINTERVALS = 1000;

struct iparams_outer_ell_stoch
{
   double mu;
   double p_red[DIMRED];
   double h_red[DIMRED];
};

// name OF FUNCTION: re_integrand_B_stoch
//
// PURPOSE
// =======
// Given the periodic point $p$ and the homoclinic point $h$ of the circular
// RTBP, an integration time $s$ in the homoclinic and periodic trajectory, 
// this function computes the integrand of the function $B^+$ or $C^+$
// (depending if this function is called from re_B_stoch or re_C_stoch),
// \[ f(\gamma_h(s)) e^{it(s)} - f(\gamma_p(s)) e^{i(t(s)+\omega_+^j)} \]
// \[ f(\gamma_h(s)) e^{it(s)} - f(\gamma_p(s)) e^{i(t(s)+\omega_-^j)} \]
//
// PARAMETERS
// ==========
// s
//    integration time in the homoclinic trajectory.
// mu
//    mass parameter for the RTBP
// p
//    periodic point, 6 coordinates: (l_p,L_p,g_p,G_p,t_p,I_p). 
// h
//    homoclinic point, 6 coordinates: (l_h,L_h,g_h,G_h,t_h,I_h). 
// 
// RETURN VALUE
// ============
// Returns the REAL part of the integrand.
//
// NOTES
// =====
// Notice that $p$ is actually a point in the periodic trajectory $\gamma_3$
// (or $\gamma_4$), and $h$ is a point in the homoclinic trajectory
// $\gamma^f$.
//
// CALLS TO: frtbp_red_g, re_f_integrand_stoch, im_f_integrand_stoch

double re_integrand_B_stoch(double s, void *params)
{
   double mu;
   double p[DIMRED];
   double h[DIMRED];

   // auxiliary variables
   int status;
   double p2[DIM], h2[DIM];

   double t_p, t_h;     // original time

   double re_fp, im_fp; 	// f(\gamma_p(s))
   double re_fh, im_fh; 	// f(\gamma_h(s))

   double term1, term2; 	// first,second term in integrand

   mu = ((struct iparams_outer_ell_stoch *)params)->mu;

   dblcpy(p, ((struct iparams_outer_ell_stoch *)params)->p_red, DIMRED);
   dblcpy(h, ((struct iparams_outer_ell_stoch *)params)->h_red, DIMRED);

   // Compute $\Phi_s(p)$
   status = frtbp_red_g(mu,s,p);
   if(status)
   {
      fprintf(stderr, "integrand: error integrating trajectory");
      exit(EXIT_FAILURE);
   }
   t_p = p[4];

   // Compute $\Phi_s(h)$
   status = frtbp_red_g(mu,s,h);
   if(status)
   {
      fprintf(stderr, "integrand: error integrating trajectory");
      exit(EXIT_FAILURE);
   }
   t_h = h[4];

   // p2=\Phi_s(p)
   dblcpy(p2,p,DIM);

   // h2=\Phi_s(h)
   dblcpy(h2,h,DIM);

   // Compute first term in integrand: f(\gamma_h(s)) e^{it_h} (real part).
   re_fh = re_f_integrand_stoch(mu,h2);
   im_fh = im_f_integrand_stoch(mu,h2);
   term1 = -(re_fh*sin(t_h) + im_fh*cos(t_h));

   // Compute second term in integrand: f(\gamma_p(s)) e^{i(t_p+\omega)}
   // (real part). Notice that t_p has already been shifted by \omega in
   // function re_B_stoch, so no need to do it here.
   re_fp = re_f_integrand_stoch(mu,p2);
   im_fp = im_f_integrand_stoch(mu,p2);
   term2 = -(re_fp*sin(t_p) + im_fp*cos(t_p));

   return term1 - term2;
}

double im_integrand_B_stoch(double s, void *params)
{
   double mu;
   double p[DIMRED];
   double h[DIMRED];

   // auxiliary variables
   int status;
   double p2[DIM], h2[DIM];

   double t_p, t_h;     // original time

   double re_fp, im_fp; 	// f(\gamma_p(s))
   double re_fh, im_fh; 	// f(\gamma_h(s))

   double term1, term2; 	// first,second term in integrand

   mu = ((struct iparams_outer_ell_stoch *)params)->mu;

   dblcpy(p, ((struct iparams_outer_ell_stoch *)params)->p_red, DIMRED);
   dblcpy(h, ((struct iparams_outer_ell_stoch *)params)->h_red, DIMRED);

   // Compute $\Phi_s(p)$
   status = frtbp_red_g(mu,s,p);
   if(status)
   {
      fprintf(stderr, "integrand: error integrating trajectory");
      exit(EXIT_FAILURE);
   }
   t_p = p[4];

   // Compute $\Phi_s(h)$
   status = frtbp_red_g(mu,s,h);
   if(status)
   {
      fprintf(stderr, "integrand: error integrating trajectory");
      exit(EXIT_FAILURE);
   }
   t_h = h[4];

   // p2=\Phi_s(p)
   dblcpy(p2,p,DIM);

   // h2=\Phi_s(h)
   dblcpy(h2,h,DIM);

   // Compute first term in integrand: f(\gamma_h(s)) e^{it_h} (imaginary part).
   re_fh = re_f_integrand_stoch(mu,h2);
   im_fh = im_f_integrand_stoch(mu,h2);
   term1 = re_fh*cos(t_h) - im_fh*sin(t_h);

   // Compute second term in integrand: f(\gamma_p(s)) e^{i(t_p+\omega)}
   // (imaginary part). Notice that t_p has already been shifted by
   // \omega in function re_B_stoch, so no need to do it here.
   re_fp = re_f_integrand_stoch(mu,p2);
   im_fp = im_f_integrand_stoch(mu,p2);
   term2 = re_fp*cos(t_p) - im_fp*sin(t_p);

   return term1 - term2;
}

// name OF FUNCTION: re_B_stoch
// CREDIT: 
//
// PURPOSE
// =======
// Given an energy level $H$, compute the integral $B^+$.
//
// \[ B^+(I) = i\int_0^{2N\pi} 
//    f(\gamma_h(s)) e^{it(s)} - f(\gamma_p(s)) e^{i(t(s)+\omega_+^j)} ds, \]
//
// where 
//    - $\gamma_h$ is the homoclinic trajectory of the reduced flow that
//    starts at the primary homoclinic point z. 
//    - $\gamma_p$ is the periodic trajectory that starts at the periodic
//    point p_3.
//
// Notice that time is reversed in the reduced flow with respect to the
// original flow.
// Equivalently, we compute it as
//
// \[ B^+(I) = i\int_0^{2N\pi}
//    f(\gamma_h(s)) e^{it(s)} - f(\gamma_p(s)) e^{i(t(s)+\omega_+^j)} ds,
// \]
//
// where 
//    - \gamma_h(s) = \Phi_{-(14M\pi-s)}{z_u,t_0+t_f}. 
//    - \gamma_p(s) = \Phi_s{l_p,L_p,0,G_p,t_0+\omega_+^j}
//
// Here, t_f is the evolution of the time component t as we reach z from z_u,
// with oposite sign.
//
// Notice that P refers to the Poincare map of the original flow, so z^u is
// actually in the STABLE manifold of the reduced flow.
// Equivalently, we split the integral from 0 to 2N\pi into N parts of size
// 2\pi:
//
// \[ B^+(I) = 
//    \Sum_{i=0,N-1}
//       i\int_0^{2\pi}
//          f(\gamma_h(s)) e^{it(s)} 
//          - f(\gamma_p(s)) e^{i(t(s)+\omega_+^j)} ds,
// \]
//
// where 
//    - \gamma_h(s) = \Phi_{-2M\pi+2i\pi+s}{z_u,t_0+t_f}
//    - \gamma_p(s) = \Phi_{2i\pi+s}{l_p,L_p,0,G_p,t_0+\omega_+^j}
//
// Since the periodic orbit is unstable, we escape from it after a long
// integration time. Thus, it is important to exploit the fact that
// \Phi_{2\pi}(l_p,L_p,0,G_p) = (l_p,L_p,0,G_p).
//
// This is computed using numerical integration.
//
// PARAMETERS
// ==========
// mu
//    mass parameter for the RTBP
// p
//    p=(l_p,L_p,g_p=0,G_p), periodic point of period 1, on the section g=0.
// zu
//    z_u=(l_h,L_h,g_h=0,G_h), homoclinic point z_u, on the section g=0.
// omega
//    $\omega_+^j$
// res
//    On return of this function, res contains the real part of $B^+(I)$.
// M
//    Number of poincare iterates to reach z from z_u
// N
//    Upper integration limit
// 
// RETURN VALUE
// ============
// Returns a non-zero error code to indicate an error and 0 to indicate
// success.
//
// NOTES
// =====
// 
// CALLS TO: re_integrand_B_stoch

int re_B_stoch(double mu, double p[DIM], double zu[DIM], double omega, double *res,
      int M, int N)
{
   double result, error;
   gsl_integration_workspace * w
               = gsl_integration_workspace_alloc (NINTERVALS);

   struct iparams_outer_ell_stoch params;

   int i;		// integration interval
   double result_i;	// intermediate result

   double tf;   // final time to reach z from z_u

   // auxiliary variables
   double pi_red[DIMRED]; 	// point pi = P^{-i}(p)
   double xi_red[DIMRED]; 	// point xi = P^{M-i}(z^u)
   double t;
   int status;

   params.mu = mu;
   
   gsl_function F;
   F.function = &re_integrand_B_stoch;
   F.params = &params;

   // Compute the final time t_f
   xi_red[0] = zu[0];
   xi_red[1] = zu[1];
   xi_red[2] = zu[2];
   xi_red[3] = zu[3];
   xi_red[4] = 0;              // t_0
   xi_red[5] = 0;              // I_0
   frtbp_red_g(mu, -2*M*M_PI, xi_red);
   tf = -xi_red[4];

   dblcpy(pi_red,p,DIM);
   // Initialize t,I components.
   pi_red[4] = omega;  // t_0+\omega_+^f
   pi_red[5] = 0;      // I_0

   result = 0;
   for(i=0; i<N; i++)
   {
      // homoclinic point z_u
       dblcpy(xi_red,zu,DIM);
      xi_red[4] = tf;	// t_0+t_f
      xi_red[5] = 0;	// I_0
      status = frtbp_red_g(mu,-(M-i)*2*M_PI,xi_red);
      if(status)
      {
	 fprintf(stderr, "re_B_stoch: error integrating trajectory");
	 exit(EXIT_FAILURE);
      }

      dblcpy(params.p_red, pi_red, DIMRED);
      dblcpy(params.h_red, xi_red, DIMRED);

      // Integrate integrand function by parts. 
      // Previously, we used 2M parts of size \pi. Now we use M parts of
      // size 2pi.
      // We request a absolute error of 0 and a relative error RELERROR.
      gsl_integration_qags (&F, 0, 2*M_PI, 0, RELERROR, NINTERVALS, w,
	    &result_i, &error);
      fprintf (stderr, "estimated error = % .3le\n", error);

      // periodic point p_3

      // It is important to exploit the fact that
      // \Phi_{2\pi}(l_p,L_p,0,G_p) = (l_p,L_p,0,G_p).
      dblcpy(pi_red, p, DIM);
      status = frtbp_red_g(mu,2*M_PI,pi_red);
      if(status)
      {
	 fprintf(stderr, "integrand: error integrating trajectory");
	 exit(EXIT_FAILURE);
      }

      result += result_i;
   }
   gsl_integration_workspace_free (w);

   *res = result;		// real(B^+)
   return 0;
}

int im_B_stoch(double mu, double p[DIM], double zu[DIM], double omega, double *res,
      int M, int N)
{
   double result, error;
   gsl_integration_workspace * w
               = gsl_integration_workspace_alloc (NINTERVALS);

   struct iparams_outer_ell_stoch params;

   int i;		// integration interval
   double result_i;	// intermediate result

   double tf;   // final time to reach z from z_u

   // auxiliary variables
   double pi_red[DIMRED]; 	// point pi = P^{-i}(p)
   double xi_red[DIMRED]; 	// point xi = P^{M-i}(z^u)
   double t;
   int status;

   params.mu = mu;
   
   gsl_function F;
   F.function = &im_integrand_B_stoch;
   F.params = &params;

   // Compute the final time t_f
   xi_red[0] = zu[0];
   xi_red[1] = zu[1];
   xi_red[2] = zu[2];
   xi_red[3] = zu[3];
   xi_red[4] = 0;              // t_0
   xi_red[5] = 0;              // I_0
   frtbp_red_g(mu, -2*M*M_PI, xi_red);
   tf = -xi_red[4];

   dblcpy(pi_red,p,DIM);
   // Initialize t,I components.
   pi_red[4] = omega;  // t_0+\omega_0^+
   pi_red[5] = 0;              // I_0

   result = 0;
   for(i=0; i<N; i++)
   {
      // homoclinic point z_u
       dblcpy(xi_red,zu,DIM);
      xi_red[4] = tf;	// t_0+t_f
      xi_red[5] = 0;	// I_0
      status = frtbp_red_g(mu,-(M-i)*2*M_PI,xi_red);
      if(status)
      {
	 fprintf(stderr, "re_B_stoch: error integrating trajectory");
	 exit(EXIT_FAILURE);
      }

      dblcpy(params.p_red, pi_red, DIMRED);
      dblcpy(params.h_red, xi_red, DIMRED);

      // Integrate integrand function by parts. 
      // Previously, we used 2M parts of size \pi. Now we use M parts of
      // size 2pi.
      // We request a absolute error of 0 and a relative error RELERROR.
      gsl_integration_qags (&F, 0, 2*M_PI, 0, RELERROR, NINTERVALS, w,
	    &result_i, &error);
      fprintf (stderr, "estimated error = % .3le\n", error);

      // periodic point p_3

      // It is important to exploit the fact that
      // \Phi_{2\pi}(l_p,L_p,0,G_p) = (l_p,L_p,0,G_p).
      dblcpy(pi_red,p,DIM);
      status = frtbp_red_g(mu,2*M_PI,pi_red);
      if(status)
      {
	 fprintf(stderr, "integrand: error integrating trajectory");
	 exit(EXIT_FAILURE);
      }

      result += result_i;
   }
   gsl_integration_workspace_free (w);

   *res = result;		// imaginary(B^+)
   return 0;
}

// name OF FUNCTION: re_C_stoch
// CREDIT: 
//
// PURPOSE
// =======
// Given an energy level $H$, compute the integral $C^+$.
//
// \[ C^+(I) = i\int_0^{-2N\pi} 
//    f(\gamma_h(s)) e^{it(s)} - f(\gamma_p(s)) e^{i(t(s)+\omega_-^j)} ds,
// \]
//
// where 
//    - $\gamma_h$ is the homoclinic trajectory of the reduced flow that
//    starts at the primary homoclinic point z. 
//    - $\gamma_p$ is the periodic trajectory that starts at the periodic
//    point p_4.
//
// Notice that time is reversed in the reduced flow with respect to the
// original flow.
// Equivalently, we compute it as
//
// \[ C^+(I) = i\int_0^{-2N\pi}
//    f(\gamma_h(s)) e^{it(s)} - f(\gamma_p(s)) e^{i(t(s)+\omega_-^j)} ds,
// \]
//
// where 
//    - \gamma_h(s) = \Phi_{2M\pi+s}{z_u,t_0+t_f}. 
//    - \gamma_p(s) = \Phi_s{l_p,L_p,0,G_p,t_0+\omega_0^-}
//
// Here, t_f is the evolution of the time component t as we reach z from z_s,
// with oposite sign.
//
// Notice that P refers to the Poincare map of the original flow, so z^s is
// actually in the UNSTABLE manifold of the reduced flow.
// Equivalently, we split the integral from 0 to -2N\pi into N parts of size
// 2\pi:
//
// \[ B^+(I) = 
//    \Sum_{i=0,N-1}
//       i\int_{0}^{-2\pi}
//          f(\gamma_h(s)) e^{it(s)} 
//          - f(\gamma_p(s)) e^{i(t(s)+\omega_-^j)} ds,
// \]
//
// where 
//    - \gamma_h(s) = \Phi_{2M\pi-2i\pi+s}{z_s,t_0+t_f}
//    - \gamma_p(s) = \Phi_{-2i\pi+s}{l_p,L_p,0,G_p,t_0+\omega_-^j}
//
// Since the periodic orbit is unstable, we escape from it after a long
// integration time. Thus, it is important to exploit the fact that
// \Phi_{-2\pi}(l_p,L_p,0,G_p) = (l_p,L_p,0,G_p).
//
// This is computed using numerical integration.
//
// PARAMETERS
// ==========
// mu
//    mass parameter for the RTBP
// p
//    p=(l_p,L_p,g_p=0,G_p), periodic point of period 1, on the section g=0.
// zu
//    z_s=(l_h,L_h,g_h=0,G_h), homoclinic point z_s, on the section g=0.
// omega
//    $\omega_-^j$
// res
//    On return of this function, res contains the real part of $C^+(I)$.
// M
//    Number of poincare iterates to reach z from z_s
// N
//    Upper integration limit
// 
// RETURN VALUE
// ============
// Returns a non-zero error code to indicate an error and 0 to indicate
// success.
//
// NOTES
// =====
// 
// CALLS TO: re_integrand_B_stoch

int re_C_stoch(double mu, double p[DIM], double zs[DIM], double omega, double *res,
      int M, int N)
{
   double result, error;
   gsl_integration_workspace * w
               = gsl_integration_workspace_alloc (NINTERVALS);

   struct iparams_outer_ell_stoch params;

   int i;		// integration interval
   double result_i;	// intermediate result

   double tf;   // final time to reach z from z_s

   // auxiliary variables
   double pi_red[DIMRED]; 	// point pi = P^{i}(p)
   double xi_red[DIMRED]; 	// point xi = P^{-(M-i)}(z^s)
   double t;
   int status;

   params.mu = mu;
   
   gsl_function F;
   F.function = &re_integrand_B_stoch;
   F.params = &params;

   // Compute the final time t_f
   xi_red[0] = zs[0];
   xi_red[1] = zs[1];
   xi_red[2] = zs[2];
   xi_red[3] = zs[3];
   xi_red[4] = 0;              // t_0
   xi_red[5] = 0;              // I_0
   frtbp_red_g(mu, 2*M*M_PI, xi_red);
   tf = -xi_red[4];

   dblcpy(pi_red,p,DIM);
   // Initialize t,I components.
   pi_red[4] = omega;  // t_0+\omega_-^j
   pi_red[5] = 0;              // I_0

   result = 0;
   for(i=0; i<N; i++)
   {
      // homoclinic point z_s
       dblcpy(xi_red,zs,DIM);
      xi_red[4] = tf;	// t_0+t_f
      xi_red[5] = 0;	// I_0
      status = frtbp_red_g(mu,(M-i)*2*M_PI,xi_red);
      if(status)
      {
	 fprintf(stderr, "re_C_stoch: error integrating trajectory");
	 exit(EXIT_FAILURE);
      }

      dblcpy(params.p_red, pi_red, DIMRED);
      dblcpy(params.h_red, xi_red, DIMRED);

      // Integrate integrand function by parts. 
      // Previously, we used 2M parts of size \pi. Now we use M parts of
      // size 2pi.
      // We request a absolute error of 0 and a relative error RELERROR.
      gsl_integration_qags (&F, 0, -2*M_PI, 0, RELERROR, NINTERVALS, w,
	    &result_i, &error);
      fprintf (stderr, "estimated error = % .3le\n", error);

      // periodic point p_4

      // It is important to exploit the fact that
      // \Phi_{-2\pi}(l_p,L_p,0,G_p) = (l_p,L_p,0,G_p).
      dblcpy(pi_red,p,DIM);
      status = frtbp_red_g(mu,-2*M_PI,pi_red);
      if(status)
      {
	 fprintf(stderr, "integrand: error integrating trajectory");
	 exit(EXIT_FAILURE);
      }

      result += result_i;
   }
   gsl_integration_workspace_free (w);

   *res = result;		// real(C^+)
   return 0;
}

int im_C_stoch(double mu, double p[DIM], double zs[DIM], double omega, double *res,
      int M, int N)
{
   double result, error;
   gsl_integration_workspace * w
               = gsl_integration_workspace_alloc (NINTERVALS);

   struct iparams_outer_ell_stoch params;

   int i;		// integration interval
   double result_i;	// intermediate result

   double tf;   // final time to reach z from z_s

   // auxiliary variables
   double pi_red[DIMRED]; 	// point pi = P^{i}(p)
   double xi_red[DIMRED]; 	// point xi = P^{-(M-i)}(z^s)
   double t;
   int status;

   params.mu = mu;
   
   gsl_function F;
   F.function = &im_integrand_B_stoch;
   F.params = &params;

   // Compute the final time t_f
   xi_red[0] = zs[0];
   xi_red[1] = zs[1];
   xi_red[2] = zs[2];
   xi_red[3] = zs[3];
   xi_red[4] = 0;              // t_0
   xi_red[5] = 0;              // I_0
   frtbp_red_g(mu, 2*M*M_PI, xi_red);
   tf = -xi_red[4];

   dblcpy(pi_red,p,DIM);
   // Initialize t,I components.
   pi_red[4] = omega;  // t_0+\omega_0^-
   pi_red[5] = 0;      // I_0

   result = 0;
   for(i=0; i<N; i++)
   {
      // homoclinic point z_s
       dblcpy(xi_red,zs,DIM);
      xi_red[4] = tf;	// t_0+t_f
      xi_red[5] = 0;	// I_0
      status = frtbp_red_g(mu,(M-i)*2*M_PI,xi_red);
      if(status)
      {
	 fprintf(stderr, "re_C_stoch: error integrating trajectory");
	 exit(EXIT_FAILURE);
      }

      dblcpy(params.p_red, pi_red, DIMRED);
      dblcpy(params.h_red, xi_red, DIMRED);

      // Integrate integrand function by parts. 
      // Previously, we used 2M parts of size \pi. Now we use M parts of
      // size 2pi.
      // We request a absolute error of 0 and a relative error RELERROR.
      gsl_integration_qags (&F, 0, -2*M_PI, 0, RELERROR, NINTERVALS, w,
	    &result_i, &error);
      fprintf (stderr, "estimated error = % .3le\n", error);

      // periodic point p_4

      // It is important to exploit the fact that
      // \Phi_{-2\pi}(l_p,L_p,0,G_p) = (l_p,L_p,0,G_p).
      dblcpy(pi_red,p,DIM);
      status = frtbp_red_g(mu,-2*M_PI,pi_red);
      if(status)
      {
	 fprintf(stderr, "integrand: error integrating trajectory");
	 exit(EXIT_FAILURE);
      }

      result += result_i;
   }
   gsl_integration_workspace_free (w);

   *res = result;		// imaginary(C^+)
   return 0;
}

