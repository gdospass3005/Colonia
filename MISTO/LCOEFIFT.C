#include<stdio.h>
#include<math.h>
#include"prisft.h"

#define EFSCAN(j)  fscanf(coef,"%lf",&j)
#define EFSCANd(j)  fscanf(coef,"%d",&j)

FILE *coef;

void load_coef_ift(char *Arq)
{
  extern long MAX_ITER;
  extern double rho,zzero,delta,HMAX;

 if ((coef = fopen(Arq,"r")) == NULL)
  {
   fprintf(stderr,"Erro na abertura do arquivo %s",Arq);
   exit(-1);
  }

 EFSCAN(rho); 
 EFSCAN(zzero);
 EFSCAN(delta);
 EFSCANd(MAX_ITER);
 EFSCAN(HMAX);

 fclose(coef);

}
