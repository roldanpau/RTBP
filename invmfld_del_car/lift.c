/*! \file
  \brief Lift points in linear segment from \R^2 to \R^4
  */

#include <stdio.h>
#include <rtbp.h>       // DIM
#include <section.h>
#include <hinv.h>

int lift(double mu, section_t sec, double H, int n, const double *l, 
        double *l4)
{
    // Auxiliary variables
    int i, status;
    const double *p; 
    double *p4;

    for(i=0; i<n; i++)
    {
        p=l+2*i;
        p4=l4+DIM*i;

        p4[0]=p[0]; // x
        p4[1]=0;    // y
        p4[2]=p[1]; // p_x
        status=hinv(mu,sec,H,p4);
       if(status)
       {
          fprintf(stderr, "lift: error lifting point\n");
          return(1);
       }
    }
    return(0);
}
