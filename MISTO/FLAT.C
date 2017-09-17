#include<stdio.h>
#include<math.h>
#include"prisft.h"

void flat(double **G, double **H)
/* Primeira estimativa da topografia --- Profundidades de Platos de Bouguer.
   Anomalias em mGal, Profundidades em metros. */
{
 extern double rho;
 extern long td_x,td_y;
 
 long i,j,k,l;
 double pi=(double)3.14159265358979;
 double gamma=(double)6.67259e-3;

 for(i=1;i<(td_x+1);i++)
  for(j=1;j<(td_y+1);j++)
    H[i][j] = G[i][j]/( (double)2. * pi * gamma * rho );

}
