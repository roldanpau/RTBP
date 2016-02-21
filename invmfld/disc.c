// name OF FUNCTION: disc
// CREDIT: 
//
// DESCRIPTION
// ===========
// Discretize the linear segment between $p_0$ and $p_1$ into NPOINTS.
//
// PARAMETERS
// ==========
// p0,p1
//    endpoints of linear segment
// n
//    number of points in the discretization
// l
//    On return of the this function, it holds the discretization of the
//    linear segment, i.e. an array of n 2D points.
// 
// RETURN VALUE
// ============
// Returns a non-zero error code to indicate an error and 0 to indicate
// success.
//
// NOTES
// =====
// n must be >1.
//
// Caller must make sure that enough space has been allocated for l (2*n
// doubles)
//
// The points are equally spaced. 
// The first point coincides with the left endpoint of the linear segment.
// The last point coincides with the right endpoint of the linear segment.

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
