#include <stdio.h>
#include "taylor2dv.h"

double mu;

int main(void)
{
  double t,tf,x[20],h;
  int nt,j;
  FILE *f;
  mu=0.01;
  t=0.e0;
  x[0]=-0.45;
  x[1]= 0.80;
  x[2]=-0.80;
  x[3]=-0.45;

  x[ 4]=1; x[ 5]=0; x[ 6]=0; x[ 7]=0;
  x[ 8]=0; x[ 9]=1; x[10]=0; x[11]=0;
  x[12]=0; x[13]=0; x[14]=1; x[15]=0;
  x[16]=0; x[17]=0; x[18]=0; x[19]=1;

  tf=1.e0;

  while (taylor_step_rtbp2dv(&t,x,1,1,-16,-16,&tf,&h,&nt,NULL) != 1);

  f=fopen("prova2dv.res","w");
  for (j=0; j<20; j++) fprintf(f,"%3d %20.15f\n",j,x[j]);
  fclose(f);
  return 0;
}
