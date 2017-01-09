/*! \file
    \brief Number of Poincare iterates to reach homoclinic point
    \author Pau Roldan
*/

#include <stdio.h>
#include <stdlib.h>	    // EXIT_SUCCESS, EXIT_FAILURE
#include <math.h>	    // fabs, fmax
#include <frtbp.h>      // DIM
#include <prtbp_2d.h>
#include <sec2sec1_module.h>

double dist(double x[2], double y[2])
{
    return fmax(fabs(x[0]-y[0]), fabs(x[1]-y[1]));
}

/** 
   Number of Poincare iterates to reach homoclinic point

   Consider the Restricted Three Body Problem. Assume that energy is fixed:
   $H(x,y,p_x,p_y)=\bar H$.
   Let $SEC2={y=0, p_y<0} be a Poincare section, and $P(x,p_x)$ be 
   the associated 2D Poincare map.
   Let $p(x,p_x)$ be a fixed point of the Poincare map associated to a
   periodic orbit of the flow.
   Let $W^u(p), W^s(p)$ be the unstable/stable manifold of $p$.
   Let $p_u$ be the preimage of the homoclinic point, located on the 
   unstable manifold.
   Let $z$ be the primary homoclinic point $z\in W^u(p) \cap W^s(p)$.
   This program computes the number of iterates of the Poincare map 
   (in Cartesian coordinates) needed to reach $z$ from $p_u$.
  
   OVERALL METHOD
   ==============
  
   1. Input parameters from stdin:
   
      - mass parameter 
      - energy value "H"
      - point "p_u" in the unstable segment (2D)
      - homoclinic point "z" (2D)
  
   2. Iterate the point p_u until it reaches homoclinic point z.
  
   3. Output the following data to stdout:

      - energy value "H"
      - section "SEC1" or "SEC2" where the homoclinic point lies
      - number of iterates to reach z
 */

int main( )
{
    const double MAX_DIST = 1e-4;
    const int MAX_ITERS= 200;
   double mu, H;

   double p_u[2];	// points in the unstable segment
   double z[DIM];	    // homoclinic points
   double z2[2];	    // homoclinic points

   // Number of iterations by the Poincare map in the unstable dir.
   int n;

   // auxiliary vars
   int status;
   double x[2];     // point in the manifold
   double ti;
   double vy;
   section_t sec;

   // 1. Input parameters from stdin.
   if(scanf("%le", &mu) < 1)
   {
      perror("main: error reading input");
      exit(EXIT_FAILURE);
   }

   while(scanf("%le %le %le %le %le %le %le", 
	    &H, p_u, p_u+1, z, z+1, z+2, z+3) == 7)
   {
       fprintf(stderr, "H=%e\n",H);

       z2[0] = z[0];    // x
       z2[1] = z[2];    // p_x

       vy = z[3]-z[0];  // p_y - x
       if(vy>0)
       {
           // take p_u to section SEC1
           sec2sec1(mu,H,p_u,&ti);
           sec = SEC1;
       }
       else
           sec = SEC2;

       x[0] = p_u[0];
       x[1] = p_u[1];
       
      // 2. Iterate the point x until it reaches homoclinic point z.
      n=0;
      while(dist(x,z2)>MAX_DIST && n<MAX_ITERS)
      {
          status = prtbp_2d(mu,sec,H,1,x,&ti);
          if(status)
          {
             fprintf(stderr, "main: error iterating point");
             exit(EXIT_FAILURE);
          }
          n++;
      }
      if(n==MAX_ITERS)
      {
         fprintf(stderr, "main: too many iterates!");
         //exit(EXIT_FAILURE);
      }

       // 3. Output the following data to stdout:
       // 
       //    - energy value "H"
       //    - section "SEC1" or "SEC2" where the homoclinic point lies
       //    - number of iterates to reach z
      printf("%le %s %d\n", H, (sec==SEC1? "SEC1":"SEC2"), n);
   }
   exit(EXIT_SUCCESS);
}
