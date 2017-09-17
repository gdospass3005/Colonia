#include<stdio.h>
#include<math.h>

void mudaorizf(float ***H,float z0)
{
 extern long tdat_x,tdat_y;
 long cont1,cont2;

 for(cont2=1;cont2<(tdat_y+1);cont2++)
  for(cont1=1;cont1<(tdat_x+1);cont1++)
   H[1][cont1][cont2] = z0 - H[1][cont1][cont2];
}
