#include <stdio.h>
#include "taylor3dv.h"

double mu;

int main(void)
{
  double t,tf,x[42],h;
  int nt,j;
  FILE *f;
  mu=0.01;
  t=0.e0;
  x[0]=-0.45;
  x[1]= 0.80;
  x[2]= 0.00;
  x[3]=-0.80;
  x[4]=-0.45;
  x[5]= 0.00;

  x[ 6]=1; x[ 7]=0; x[ 8]=0; x[ 9]=0; x[10]=0; x[11]=0;
  x[12]=0; x[13]=1; x[14]=0; x[15]=0; x[16]=0; x[17]=0;
  x[18]=0; x[19]=0; x[20]=1; x[21]=0; x[22]=0; x[23]=0;
  x[24]=0; x[25]=0; x[26]=0; x[27]=1; x[28]=0; x[29]=0;
  x[30]=0; x[31]=0; x[32]=0; x[33]=0; x[34]=1; x[35]=0;
  x[36]=0; x[37]=0; x[38]=0; x[39]=0; x[40]=0; x[41]=1;

  tf=1.e0;

  while (taylor_step_rtbp3dv(&t,x,1,1,-16,-16,&tf,&h,&nt,NULL) != 1);

  f=fopen("prova3dv.res","w");
  for (j=0; j<42; j++) fprintf(f,"%3d %20.15f\n",j,x[j]);
  fclose(f);
  return 0;
}
