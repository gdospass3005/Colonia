#include<stdio.h>
#include<math.h>
#include"prisft.h"

void mudaoriz(double **H,double z0)
{
 extern long td_x,td_y;
 long cont1,cont2;

 for(cont2=1;cont2<(td_y+1);cont2++)
  for(cont1=1;cont1<(td_x+1);cont1++)
   H[cont1][cont2] = z0 - H[cont1][cont2];
}
