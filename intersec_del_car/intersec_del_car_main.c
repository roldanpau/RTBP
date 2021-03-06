/*! \file
    \brief Intersection of invariant manifolds: main program.
    \author Pau Roldan
*/

#include <stdio.h>
#include <stdlib.h>	// EXIT_SUCCESS, EXIT_FAILURE
#include <string.h>	// strcmp
#include <rtbp.h>	// DIM

#include <section.h>    // section_t, branch_t
#include <errmfld.h>    // h_opt

#include <prtbp_nl_2d_module.h>

#include "intersec_del_car.h"

#include <utils_module.h>	// dblprint

/**
 * Intersection of Invariant Manifolds.
 *
  \param[in] mu         mass parameter for the RTBP
  \param[in] sec        Poincare section: sec={SEC1,SEC2}
  \param[in] br         branch type: br={LEFT, RIGHT}

  \param[in] stable     stable flag specifies wheather we want to compute the
  unstable (=0) or stable (=1) manifold

  \param[in] H          energy value

  \param[in] n
  number of desired iterations by the Poincare map in the unstable direction.

  \param[in] p_u[2]
  Point in the unstable segment (Cartesian coords)

  \param[out] t
  On exit, it contains integration time to reach z from p_u.
  
  \param[out] z_del[DIM]	
  On exit, it contains homoclinic point z = P^n(p_u) in Delaunay.

  \param[out] z_car[DIM]	
  On exit, it contains homoclinic point z = P^n(p_u) in Cartesian.

  \param[out] z_u[DIM]
  Point in local unstable manifold of the appropriate pendulum (Delaunay
  coords). This will be needed in outer_circ.
  
  \param[out] z_u_car[DIM]
  Point in local unstable manifold of the appropriate pendulum (Cartesian
  coords). This will be needed in outer_circ.
  
  \returns 
  a non-zero error code to indicate an error and 0 to indicate
  success.

 * OVERALL METHOD
 * ==============
 *
 * 1. Input parameters from stdin:
 * 
 *    - mass parameter 
 *    - section type "sec": sec={SEC1,SEC2}
 *    - branch type "br": br={LEFT, RIGHT}
 *    - stable flag
 *    - axis line "l"
 *
 * For each energy value, do:
 *
 * 2.1 Input data:
 *    [1] energy value "H"
 *    [2-3] fixed point "p"
 *    [4-5] linear unstable direction "v"
 *    [5] unstable eigenvalue "lambda"
 *    [6] number of iterations "n" in the unstable dir by the Poincare map to
 *     reach the intersection point
 *    [7-8] "h1, h2" small increment in the direction of v (initial guess for
 *    bisection).
 *
 * 2.2. Find a root of the distance function, i.e. an intersection point of the
 * manifold with the axis line.
 *
 * 2.3. Output the following data to stdout:
 *    [1] energy level H
 *    [2-3] point p_u,
 *    [4-4] t
 *      integration time to reach z from p_u.
 *    [5-8] z_del[DIM]	
 *      homoclinic point z = P^n(p_u) in Delaunay.
 *    [9-12] z_car[DIM]	
 *      homoclinic point z = P^n(p_u) in Cartesian.
 *    [13-16] z_u[DIM]
 *      Point in local unstable manifold of the appropriate pendulum
 *      (Delaunay coords). This will be needed in outer_circ.
 *    [17-20] z_u_car[DIM]
 *      Point in local unstable manifold of the appropriate pendulum
 *      (Cartesian coords). This will be needed in outer_circ.
 */

int main( )
{
   double mu, H;
   section_t sec;   // Poincare section
   branch_t br;     // branch type

   double p[2];		// fixed point
   double v[2];		// eigenvector
   double lambda;	// eigenvalue

   // number of desired iterations by the Poincare map
   int n;		// in the unstable dir

   // Small increment in the direction of v (initial guess for Newton).
   // This is furnished by program "approxint".
   double h1, h2;

   // "stable" flag specifies wheather we want to compute the unstable (=0)
   // or stable (=1) manifold
   int stable;

   double l;		// axis line p_x = l

   double p0[2];        // p0 = p+hv
   double p1[2];        // p1 = P(p0)

   double h;		// root of distance function

   double p_u[2];	// point in the unstable segment
   double z[2];		// homoclinic point

   double z_car[DIM];		// homoclinic point in Cartesian
   double z_del[DIM];		// homoclinic point in Delaunay
   double z_u[DIM];		    // point in local unstable manifold (Delaunay)
   double z_u_car[DIM];		// point in local unstable manifold (Cartesian)

   double t;		// integration time to reach z from p_u/p_s

   // auxiliary vars
   char section_str[10];    // holds input string "SEC1", "SEC2" etc
   char branch_str[10];    // holds input string "LEFT", "RIGHT" etc
   int status;

   // 1. Input parameters from stdin.
   if(scanf("%le %s %s %d %le", &mu, section_str, branch_str, &stable, &l) < 5)
   {
      perror("main: error reading input");
      exit(EXIT_FAILURE);
   }

   if (strcmp(section_str,"SEC1") == 0)
      sec = SEC1;
   else if (strcmp(section_str,"SEC2") == 0)
      sec = SEC2;
   else if (strcmp(section_str,"SECg") == 0)
      sec = SECg;
   else if (strcmp(section_str,"SECg2") == 0)
      sec = SECg2;
   else
   {
      perror("main: error reading section string");
      exit(EXIT_FAILURE);
   }

   if (strcmp(branch_str,"RIGHT") == 0)
      br = RIGHT;
   else if (strcmp(branch_str,"LEFT") == 0)
      br = LEFT;
   else
   {
      perror("main: error reading branch string");
      exit(EXIT_FAILURE);
   }

   while(scanf("%le %le %le %le %le %le %d %le %le", 
	    &H, p, p+1, v, v+1, &lambda, &n, &h1, &h2) == 9)
   {
      fprintf(stderr, "Processing energy level %le...\n", H);

	  h=1.e-2;	// upper bound on linear displacement

      // Swap branches, to have RIGHT branch = upper branch in Delaunay, and
      // LEFT branch = lower branch in Delaunay.
      // This is only necessary in the stable case.
      if(stable)
          h=-h;

      // By default we work with the RIGHT branch of the manifolds.
      // To work with the LEFT branch of the manifolds instead,
      // we just take the negative of the displacement h in the linear
      // approximation.
      if(br==LEFT) h=-h;

      // Compute optimal displacement h such that the estimate error commited
      // in the linear approximation of the manifold is smallest.
      h = h_opt(mu,SEC2,H,4,p,v,lambda,stable,h);

   // Compute $p_0$
   p0[0] = p[0] + h*v[0];
   p0[1] = p[1] + h*v[1];

   // Compute $p_1$
   p1[0] = p0[0];
   p1[1] = p0[1];
   status=prtbp_nl_2d(mu,SEC2,H,4,p1,&t);        // $p_1 = P(p_0)$
   if(status)
   {
      fprintf(stderr, "approxint_del_car: error computing Poincare map\n");
      return(1);
   }

      // 2. Find a root of the distance function, i.e. an intersection point of
      // the manifolds
      if(!stable)
      {
          status = intersec_del_car_unst(mu, sec, br, H, p0, p1, lambda, n, 
                  h1, h2, l, &h, p_u, &t, z_del, z_car, z_u, z_u_car);
          if(status)
          {
              fprintf(stderr, "main: error computing intersection point\n");
              exit(EXIT_FAILURE);
          }
      }
      else
      {
          status = intersec_del_car_st(mu, sec, br, H, p, v, lambda, n, 
                  h1, h2, l, &h, p_u, &t, z_del, z_car, z_u, z_u_car);
          if(status)
          {
              fprintf(stderr, "main: error computing intersection point\n");
              exit(EXIT_FAILURE);
          }
      }

      // 3. Output the following data to stdout:
      //    - energy level H
      //    - point p_u/p_s
      //    - integration time t to reach the intersection point z, 
      //    - intersection point z = P(p_u).
      //    - point in local unstable manifold z_u
      //    - point in local unstable manifold z_u_car
      printf("%.15le ", H);
      dblprint(p_u,2);
      printf("%.15le ", t);
      dblprint(z_del,DIM);
      dblprint(z_car,DIM);
      dblprint(z_u,DIM);
      dblprint(z_u_car,DIM);
      printf("\n");
      fflush(NULL);

      fprintf(stderr, "Done!\n");
      fflush(stderr);
   }
   exit(EXIT_SUCCESS);
}
