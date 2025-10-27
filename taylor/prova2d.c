#include <stdio.h>
#include "taylor2d.h"

double mu;

int main(void)
{
  double t,tf,x[4],h;
  int nt,j;
  FILE *f;
  mu=0.01;
  t=0.e0;
  x[0]=-0.45;
  x[1]= 0.80;
  x[2]=-0.80;
  x[3]=-0.45;

  tf=1.e0;

  while (taylor_step_rtbp2d(&t,x,1,1,-16,-16,&tf,&h,&nt,NULL) != 1);

  f=fopen("prova2d.res","w");
  for (j=0; j<4; j++) fprintf(f,"%3d %20.15f\n",j,x[j]);
  fclose(f);
  return 0;
}
