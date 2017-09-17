#include<stdio.h>
#include<math.h>
#include"prisver.h"

#define EFSCANd(j)  fscanf(entrada,"%D",&j)
#define EFSCAN(j)  fscanf(entrada,"%lf",&j)
#define EFSCANs(j) fscanf(entrada,"%s",j)

FILE *entrada;
char cod[10];

void load_g(char *Arq)
{
  extern double **g,interv_x,interv_y,difinterv,xmin,xmax,ymin,ymax,gmin,gmax;
  extern long td_x,td_y;
  long cont1,cont2;

 if ((entrada = fopen(Arq,"r")) == NULL)
  {
   fprintf(stderr,"Erro na abertura do arquivo %s",Arq);
   exit(-1);
  }

 EFSCANs(cod);
 EFSCANd(td_x); EFSCANd(td_y);
 EFSCAN(xmin);
 EFSCAN(xmax);
 EFSCAN(ymin);
 EFSCAN(ymax);
 EFSCAN(gmin);
 EFSCAN(gmax);

 interv_x = (xmax - xmin)/(double)(td_x - 1);
 interv_y = (ymax - ymin)/(double)(td_y - 1);

 difinterv = interv_x - interv_y;
 difinterv = difinterv * difinterv;
 if(difinterv > (double)0.00000001)
  {
   printf("ATENCAO: Grid deve ter espacamento igual em x e y.\n");
   printf("Saindo do sistema.\n");
   exit(-1);
  }

 g = dmatrix((long)1,td_x,(long)1,td_y);

 for(cont2=1;cont2<(td_y+1);cont2++)
   for(cont1=1;cont1<(td_x+1);cont1++)
     EFSCAN(g[cont1][cont2]);

 fclose(entrada);
}
