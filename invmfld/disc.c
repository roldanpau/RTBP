/*! \file
  \brief Discretize a linear segment into a given number of points
  */

int disc(double p0[2], double p1[2], int n, double *l)
{
   double delta[2];
   double v[2];			// v = p1-p0

   // auxiliary variables
   int i;

   // Discretize the linear segment between $p_0$ and $p_1$ into NPOINTS
   v[0] = p1[0]-p0[0];
   v[1] = p1[1]-p0[1];
   delta[0] = v[0]/(n-1);
   delta[1] = v[1]/(n-1);
   for(i=0; i<n-1; i++)
   {
      l[2*i] = p0[0]+i*delta[0];
      l[2*i+1] = p0[1]+i*delta[1];
   };
   // Last point
   l[2*i] = p1[0];
   l[2*i+1] = p1[1];
   return(0);
}
