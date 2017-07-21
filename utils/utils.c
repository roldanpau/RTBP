#include "utils_module.h"

int main()
{
   double x[] = {1, 2, 3, 4};
   double y[4];

   dblprint(x,4);

   dblcpy(y,x,4);
   dblprint(y,4);
}
