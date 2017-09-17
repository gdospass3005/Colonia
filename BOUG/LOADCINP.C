#include<stdio.h>
#include<math.h>
#include"prisboug.h"

#define EFSCAN(j)  fscanf(coef,"%lf",&j)
#define EFSCANd(j)  fscanf(coef,"%d",&j)

FILE *coef;

void load_coef(char *Arq)
{
  extern long MAX_ITER,ITER1,ITER2;
  extern double a0,a1,a2,HMAX;

 if ((coef = fopen(Arq,"r")) == NULL)
  {
   fprintf(stderr,"Erro na abertura do arquivo %s",Arq);
   exit(-1);
  }

 EFSCAN(a0);
 EFSCAN(a1);
 EFSCAN(a2);
 EFSCANd(MAX_ITER);
 EFSCANd(ITER1); 
 EFSCANd(ITER2);
 EFSCAN(HMAX);

 fclose(coef);

}
