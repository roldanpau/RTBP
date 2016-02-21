// ===================================================================
// Estimate error commited in the linear approximation of the manifold
// ===================================================================
// FILE:          $RCSfile: errmflddel.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2012-06-06 10:20:10 $
//
// PURPOSE
// =======
//
// OVERALL METHOD
// ==============

// NOTES
// =====
// If the flag "stable" specifies the unstable manifold (0), we iterate the
// forward Poincare map $P$.
// If the flag "stable" specifies the stable manifold (1), we iterate the
// backward Poincare map $P^{-1}$.
//
// We assume that manifolds are on the SEC2 section {l=pi}.

#include <stdio.h>
#include <stdlib.h>	// EXIT_SUCCESS
#include <math.h>	// sqrt
#include <prtbpdel_2d.h>	// prtbp_del_2d, prtbp_del_2d_inv

double err_mfld_del(double mu, double H, int k, double p[2], double v[2], 
      double lambda, int stable, double h)
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
      status=prtbp_del_2d(mu,SEC2,H,k,p1,&ti); 	// $p_1 = P(p_0)$
   else 	// stable manifold
      status=prtbp_del_2d_inv(mu,SEC2,H,k,p1,&ti);	// $p_1 = P^{-1}(p_0)$
   if(status)
   {
      fprintf(stderr, "err_mfld_del: error computing Poincare map\n");
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

// Compute optimal h
double h_opt_del(double mu, double H, int k, double p[2], double v[2], 
      double lambda, int stable)
{
   double h;
   double err, err2;

   h=1.e-3;
   err2 = err_mfld_del(mu,H,k,p,v,lambda,stable,h);
   do
   {
      err = err2;
      h = h/10.0;
      err2=err_mfld_del(mu,H,k,p,v,lambda,stable,h);
   } while((err/err2)>10);

   fprintf(stderr,"Optimal displacement: %e\n", h*10.0);
   fprintf(stderr,"Estimated error of manifold: %e\n", err);
   return h*10.0;
}
