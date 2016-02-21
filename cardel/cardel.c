// Headers
#include <stdio.h>	// fprintf
#include <math.h>	// sqrt
#include <rtbp.h>	// DIM

void cardel(double X[DIM], double Y[DIM])
{
   double x=X[0];
   double y=X[1];
   double px=X[2];
   double py=X[3];

   double l,L,g,G;

   // auxiliary variables
   double r = sqrt(x*x+y*y);	// radius (polar coords)

   // Note: the polar angle phi will jump from time to time 
   // e.g. when (x,y) changes from the 2nd quadrant to the 3rd
   // then phi jumps from pi to -pi. 
   // But this is inevitable, since we don't want the polar angle 
   // to grow forever.
   double phi = atan2(y,x); //+ M_PI;	// angle (polar coords)
   double dotr = px*cos(phi)+py*sin(phi); 	// radial velocity

   double e;	// eccentricity
   double u, v;

   // auxiliary variables
   double cu;	// cos(u)

   G = -y*px + x*py;
   L = sqrt( -1.0/( px*px + py*py -2.0/r ));

   // L could also be computed like this, where J=energy.
   // To avoid using an extra parameter (J), we prefer the previous
   // expression.
   // L = sqrt( -0.5/( J+G-mDH ) );

   e = sqrt(1.0 - G*G/(L*L));

   cu = 1.0/e*(1.0-r/(L*L));

   // Make sure cos(u) is in the range [-1,1].
   // It can sometimes be outside the range due to roundoff.
   if(cu<-1)
   {
      fprintf(stderr, "warning: cos(u)=%e outside range [-1,1]\n", cu);
      cu = -1;
   }
   else if(cu>1)
   {
      fprintf(stderr, "warning: cos(u)=%e outside range [-1,1]\n", cu);
      cu = 1;
   }

   if(dotr>=0)
   {
      // we are leaving the perihelion, so eccentric anomaly $u>0$
      u = acos(cu);
   }
   else
   {
      // we are approaching the perihelion, so eccentric anomaly $u<0$
      u = -acos(cu);
   }

   l = u-e*sin(u);

   if(u==M_PI)
   {
      // tan(u/2) blows up! But in this case $v=pi$.
      v=M_PI;
   }
   else
      v = 2.0*atan( sqrt((1.0+e)/(1.0-e))*tan(u/2.0) );

   // We need to determine the sign of v in [-pi,pi]
   if(0<=u && u<=M_PI)
   {
      // The sign of v must be non-negative
      if(v<0) v += M_PI;
   }
   else if(-M_PI<=u && u<0)
   {
      // The sign of v must be negative
      if(v>0) v -= M_PI;
   }

   g = phi - v;

   // On output, we normalize the angles $l,g$ between [0,2\pi).
   if(l<0) l+=2*M_PI;
   if(g<0) g+=2*M_PI;

   Y[0]=l;
   Y[1]=L;
   Y[2]=g;
   Y[3]=G;
}
