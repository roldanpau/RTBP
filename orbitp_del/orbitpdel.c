// ======================================
// Orbit of Poincare map in Delaunay RTBP
// ======================================
// FILE:          $RCSfile: orbitpdel.c,v $
// AUTHOR:        $Author: roldan $
//                All code is my own except where credited to others.
// DATE:          $Date: 2012-05-02 10:03:50 $
//
// FUNCTIONS
// =========
//
// orbitp_del
// ----------
// Consider the Restricted Three Body Problem in Delaunay variables.
// Let ${g=0} be a Poincare section, and $P(l,L,g,G)$ be the associated
// Poincare map.
// Given a point $p=(l,L,g,G)$, this function computes the orbit ${P^k(p)}$
// for $k=1,2,\dots,n$.
//
// orbitp_del_bwd
// --------------
// Consider the Restricted Three Body Problem in Delaunay variables.
// Let ${g=0} be a Poincare section, and $P(l,L,g,G)$ be the associated
// Poincare map.
// Given a point $p=(l,L,g,G)$, this function computes the backward orbit
// ${P^{-k}(p)}$ for $k=1,2,\dots,n$.
//
// print_orbit_del
// ---------------
// Print orbit of a point to stdout

#include <stdio.h>	// fprintf
#include <stdlib.h>	// EXIT_FAILURE
#include <math.h>	// M_PI
#include <prtbpdel.h>	// prtbp_del, prtbp_del_inv

// name OF FUNCTION: orbitp_del
// CREDIT: 
//
// PURPOSE
// =======
// Consider the Restricted Three Body Problem in Delaunay variables.
// Let ${g=0} be a Poincare section, and $P(l,L,g,G)$ be the associated
// Poincare map.
// Given a point $p=(l,L,g,G)$, this function computes the orbit ${P^k(p)}$
// for $k=1,2,\dots,n$.
//
// PARAMETERS
// ==========
// mu
//    mass parameter for the RTBP
// p
//    initial point of the orbit
// n
//    Last iterate of Poincare map is $P^n(p)$. The number of points in the
//    orbit is $n$.
// orbit
//    On exit, "orbit" contains the orbit of $p$. It is an array of $n$
//    consecutive 4D points: $(P^1(p), P^2(p), P^3(p), \cdots, P^n(p))$.
// 
// RETURN VALUE
// ============
// Returns a non-zero error code to indicate an error and 0 to indicate
// success.
//
// NOTES
// =====
//
// CALLS TO: prtbp_del

int orbitp_del(double mu, double p[DIM], int n, double *orbit)
{
   int i=0;
   double ti;

   // auxiliary variables
   double l, g;
   int q;	

   while(i<n)
   {
      if(prtbp_del(mu,1,p,&ti))
      {
	 fprintf(stderr, "orbitp_del: error computing Poincare map\n");
	 return(1);
      }

      // prtbp_del already normalizes $l$ and $g$ angles to $(0,2\pi)$, no
      // need to do it here.

      orbit[DIM*i] = p[0];
      orbit[DIM*i+1] = p[1];
      orbit[DIM*i+2] = p[2];
      orbit[DIM*i+3] = p[3];
      i++;
   }
   return(0);
}

// name OF FUNCTION: orbitp_del_bwd
// CREDIT: 
//
// PURPOSE
// =======
// Consider the Restricted Three Body Problem in Delaunay variables.
// Let ${g=0} be a Poincare section, and $P(l,L,g,G)$ be the associated
// Poincare map.
// Given a point $p=(l,L,g,G)$, this function computes the backward orbit
// ${P^{-k}(p)}$ for $k=1,2,\dots,n$.
//
// PARAMETERS
// ==========
// mu
//    mass parameter for the RTBP
// p
//    initial point of the orbit
// n
//    Last iterate of Poincare map is $P^{-n}(p)$. The number of points in
//    the orbit is $n$.
// orbit
//    On exit, "orbit" contains the orbit of $p$. It is an array of $n$
//    consecutive 4D points: 
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
// CALLS TO: prtbp_del_inv

int orbitp_del_bwd(double mu, double p[DIM], int n, double *orbit)
{
   int i=0;
   double ti;

   // auxiliary variables
   double l, g;
   int q;	

   while(i<n)
   {
      if(prtbp_del_inv(mu,1,p,&ti))
      {
	 fprintf(stderr, \
	       "orbitp_del_bwd: error computing inverse Poincare map\n");
	 return(1);
      }

      // prtbp_del already normalizes $l$ and $g$ angles to $(0,2\pi)$, no
      // need to do it here.

      orbit[DIM*i] = p[0];
      orbit[DIM*i+1] = p[1];
      orbit[DIM*i+2] = p[2];
      orbit[DIM*i+3] = p[3];
      i++;
   }
   return(0);
}

// name OF FUNCTION: print_orbit_del
// CREDIT: 
// PURPOSE:
// Print orbit of a point to stdout.
// The format in which the orbit is printed is:
//    p_1
//    p_2
//    p_3
//    ...
//    p_n
// where $p_k=(l,L,g,G)$ is the $k$-th point in the orbit.
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

inline void print_orbit_del(int n, double *orbit)
{
   int i;
   //for(i=0;i<n;i++)
   // Print just the last point. This is just a dirty hack in order to plot
   // the invariant manifold of a certain poincare iterate (e.g., just the
   // inv. mflds of the 1st iterate, or just the 2nd iterate, etc.)
   for(i=n-1;i<n;i++)
   {
      if(printf("% .15le % .15le % .15le % .15le\n", \
	       orbit[DIM*i], orbit[DIM*i+1], orbit[DIM*i+2],
	       orbit[DIM*i+3])<0)
      {
	 perror("print_orbit_del: error writting output");
	 exit(EXIT_FAILURE);
      }
   }
}
