/*! \file errmfld.c
    \brief Estimate error commited in the linear approximation of the manifold.
    
    $Author: roldan $
    $Date: 2013-03-11 11:19:12 $
*/

#include <stdio.h>
#include <stdlib.h>	// EXIT_SUCCESS
#include <math.h>	// sqrt
#include <assert.h>
#include <section.h>

#include <prtbp_nl_2d.h>

/** 
  Estimate error commited in the linear approximation of the manifold.

  Let $p$ be a hyperbolic fixed point for the 2D map \f$\mathcal{P}\f$.
  Assume that \f$\lambda\f$ is the unstable eigenvalue, with \f$\lambda>1\f$
  (case stable flag=0), or \f$\lambda\f$ is the stable eigenvalue, with
  \f$\lambda<1\f$ (case stable flag=1).
  Let $v$ be the unstable (resp. stable) eigenvector for the eigenvalue
  \f$\lambda\f$. 
  Given a displacement $h$ from the fixed point along the linear eigenspace,
  we consider the point $p+hv$. The error commited in the local approximation
  of the manifold is given by
  \f[ err(h) = ||\mathcal{P}(p+hv)-p-\lambda hv|| \in O(h^2). \f]

  \param[in] mu 	mass parameter for the RTBP
  \param[in] sec        type of Poincare section (sec = SEC1 or SEC2).
  \param[in] H 		energy value
  \param[in] k 		number of iterates of the Poincare map
  \param[in] p 		fixed point $p=(x,p_x)$
  \param[in] v 		eigenvector associated to unstable (or stable) direction
  \param[in] lambda 	eigenvalue associated to unstable (or stable) direction

  \param[in] 
  stable flag indicating unstable manifold (stable=0) or stable manifold
  (stable=1).

  \param[in] h linear displacement along the manifold.

  \return Returns the error commited in the linear approximation of the
  manifold.
*/

double err_mfld(double mu, section_t sec, double H, int k, double p[2],
      double v[2], double lambda, int stable, double h)
{
   double p1[2];	// p1 = P(p+hv)
   double e[2];
   double err;

   // Auxiliary variables
   int status;
   double ti;

   // Compute $p_1$
   p1[0] = p[0] + h*v[0];
   p1[1] = p[1] + h*v[1];
   if(!stable) 	// unstable manifold
      status=prtbp_nl_2d(mu,sec,H,k,p1,&ti); 	// $p_1 = P(p_0)$
   else 	// stable manifold
      status=prtbp_nl_2d_inv(mu,sec,H,k,p1,&ti);	// $p_1 = P^{-1}(p_0)$
   if(status)
   {
      fprintf(stderr, "err_mfld: error computing Poincare map\n");
      exit(EXIT_FAILURE);
   }

   if(!stable) // unstable manifold
   {
      e[0] = p1[0]-p[0]-lambda*h*v[0];
      e[1] = p1[1]-p[1]-lambda*h*v[1];
   }
   else
   {
      e[0] = p1[0]-p[0]-1.0/lambda*h*v[0];
      e[1] = p1[1]-p[1]-1.0/lambda*h*v[1];
   }
   err = sqrt(e[0]*e[0]+e[1]*e[1]);
   return(err);
}

/** 
  Compute optimal displacement $h$ from the fixed point along the linear
  eigenspace.

  Let $p$ be a hyperbolic fixed point for the 2D map \f$\mathcal{P}\f$.
  Assume that \f$\lambda\f$ is the unstable eigenvalue, with \f$\lambda>1\f$
  (case stable flag=0), or \f$\lambda\f$ is the stable eigenvalue, with
  \f$\lambda<1\f$ (case stable flag=1).
  Let $v$ be the unstable (resp. stable) eigenvector for the eigenvalue
  \f$\lambda\f$. 
  Let $h$ be a displacemnt from the fixed point along the linear eigenspace.
  If $h$ is chosen too large, then the local error approximating the manifold
  is too large. On the other hand, if $h$ is chosen too small, then we need
  to iterate the map \f$\mathcal{P}\f$ in order to globalize the manifold, and
  this leads to accumulation of numerical errors. 

  This function computes an ``optimal'' value of $h$, according to the
  following algorithm: $h$ is iteratively decreased by a factor 0.1 until the
  local error is not reduced much.
  We also harcode a minimum threshold for the local error, $err(h)>10^{-8}$,
  because we don't want to decrease $h$ too much.

NOTE: Probably a better way to choose $h$ is to minimize the empirical
numerical error commited in computing the homoclinic point.

  \param[in] mu 	mass parameter for the RTBP
  \param[in] sec        type of Poincare section (sec = SEC1 or SEC2).
  \param[in] H 		energy value
  \param[in] k 		number of iterates of the Poincare map
  \param[in] p 		fixed point $p=(x,p_x)$
  \param[in] v 		eigenvector associated to unstable (or stable) direction
  \param[in] lambda 	eigenvalue associated to unstable (or stable) direction

  \param[in] 
  stable flag indicating unstable manifold (stable=0) or stable manifold
  (stable=1).

  \return Returns optimal linear displacement along the manifold $h$.
*/

double h_opt(double mu, section_t sec, double H, int k, double p[2], 
      double v[2], double lambda, int stable)
{
   double h;
   double err, err2;

   h=1.e-2;
   err2 = err_mfld(mu,sec,H,k,p,v,lambda,stable,h);
   do
   {
      err = err2;
      h = h/10.0;
      err2=err_mfld(mu,sec,H,k,p,v,lambda,stable,h);
   } while(err2>1.e-4 && (err/err2)>10);

   fprintf(stderr,"Optimal displacement: %e\n", h*10.0);
   fprintf(stderr,"Estimated error of manifold: %e\n", err);
   return h*10.0;
}
