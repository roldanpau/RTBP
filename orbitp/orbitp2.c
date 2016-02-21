// =============================
// Orbit of Poincare map in RTBP
// =============================
// FILE:          $RCSfile: orbitp2.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2012-05-02 09:43:04 $
//
// FUNCTIONS
// =========
//
// orbitp2
// -------
// Consider the Restricted Three Body Problem. Assume that energy is fixed:
// $H(x,y,p_x,p_y)=\bar H$.
// Let $\Sigma^- = {y=0, p_y<0}$ be a Poincare section, and $P(x,p_x)$ be the
// associated 2D Poincare map.
// Let $\tilde P = P^k$ be the iterated Poincare map.
// Given a point $p=(x,p_x)$, this function computes the orbit 
// ${\tilde P^i(p)}$ for $i=1,2,\dots,n$.
//
// orbitp2_bwd
// -----------
// Consider the Restricted Three Body Problem. Assume that energy is fixed:
// $H(x,y,p_x,p_y)=\bar H$.
// Let $\Sigma^- = {y=0}$ be a Poincare section, and $P(x,p_x)$ be the
// associated 2D Poincare map.
// Given a point $p=(x,p_x)$, this function computes the backwards orbit
// ${P^{-k}(p)}$ for $k=1,2,\dots,n$.
//
// print_orbit
// -----------
// Print orbit of a point to stdout

#include <stdio.h>	// fprintf
#include <stdlib.h>	// EXIT_FAILURE

// name OF FUNCTION: orbitp2
// CREDIT: 
//
// PURPOSE
// =======
// Consider the Restricted Three Body Problem. Assume that energy is fixed:
// $H(x,y,p_x,p_y)=\bar H$.
// Let $\Sigma^- = {y=0, p_y<0}$ be a Poincare section, and $P(x,p_x)$ be the
// associated 2D Poincare map.
// Let $\tilde P = P^k$ be the iterated Poincare map.
// Given a point $p=(x,p_x)$, this function computes the orbit 
// ${\tilde P^i(p)}$ for $i=1,2,\dots,n$.
//
// PARAMETERS
// ==========
// mu
//    mass parameter for the RTBP
// H
//    energy value
// k
//    iterated Poincare map $\tilde P = P^k$ (e.g. k=6 for the 1:7 resonance)
// p
//    initial point of the orbit
// n
//    Last iterate of Poincare map is $\tilde P^n(p)$. The number of points
//    in the orbit is $n$.
// orbit
//    On exit, "orbit" contains the orbit of $p$. It is an array of $n$
//    consecutive 2D points: $(\tilde P^1(p), \tilde P^2(p), \tilde P^3(p),
//    \cdots, \tilde P^n(p))$.
// 
// RETURN VALUE
// ============
// Returns a non-zero error code to indicate an error and 0 to indicate
// success.
//
// NOTES
// =====
//
// CALLS TO: prtbp2_2d

int orbitp2(double mu, double H, int k, double p[2], int n, double *orbit)
{
   int i=0;
   double ti;

   while(i<n)
   {
      if(prtbp2_2d(mu,H,k,p,&ti))
      {
	 fprintf(stderr, "orbitp2: error computing Poincare map\n");
	 return(1);
      }
      orbit[2*i] = p[0];
      orbit[2*i+1] = p[1];
      i++;
   }
   return(0);
}

// name OF FUNCTION: orbitp2_bwd
// CREDIT: 
//
// PURPOSE
// =======
// Consider the Restricted Three Body Problem. Assume that energy is fixed:
// $H(x,y,p_x,p_y)=\bar H$.
// Let $\Sigma^- = {y=0, p_y<0}$ be a Poincare section, and $P(x,p_x)$ be the
// associated 2D Poincare map.
// Given a point $p=(x,p_x)$, this function computes the backwards orbit
// ${P^{-k(p)}$ for $k=1,2,\dots,n$.
//
// PARAMETERS
// ==========
// mu
//    mass parameter for the RTBP
// H
//    energy value
// p
//    initial point of the orbit
// n
//    Last iterate of Poincare map is $P^{-n}(p)$. The number of points in
//    the orbit is $n$.
// orbit
//    On exit, "orbit" contains the orbit of $p$. It is an array of $n$
//    consecutive 2D points: 
//       $(P^{-1}(p), P^{-2}(p), P^{-3}(p), \cdots, P^{-n}(p))$.
// 
// RETURN VALUE
// ============
// Returns a non-zero error code to indicate an error and 0 to indicate
// success.
//
// NOTES
// =====
//
// CALLS TO: prtbp2_2d_inv

int orbitp2_bwd(double mu, double H, double p[2], int n, double *orbit)
{
   int i=0;
   double ti;

   while(i<n)
   {
      if(prtbp2_2d_inv(mu,H,6,p,&ti))
      {
	 fprintf(stderr, "orbitp2_bwd: error computing inverse Poincare map\n");
	 return(1);
      }
      orbit[2*i] = p[0];
      orbit[2*i+1] = p[1];
      i++;
   }
   return(0);
}

// name OF FUNCTION: print_orbit
// CREDIT: 
// PURPOSE:
// Print orbit of a point to stdout.
// The format in which the orbit is printed is:
//    p_1
//    p_2
//    p_3
//    ...
//    p_n
// where $p_k=(x,p_x)$ is the $k$-th point in the orbit.
//
// PARAMETERS
// ==========
// n 
//    number of points in the orbit.
// orbit 
//    array of n points
// 
// RETURN VALUE:
//
// CALLS TO: none
//
// CALLED FROM: main

inline void print_orbit(int n, double *orbit)
{
   int i;
   for(i=0;i<n;i++)
   {
      if(printf("% .15le % .15le\n", orbit[2*i], orbit[2*i+1])<0)
      {
	 perror("main: error writting output");
	 exit(EXIT_FAILURE);
      }
   }
}
