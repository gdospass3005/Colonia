#include<stdio.h>
#include<math.h>
#include"prisft.h"

#define EFSCANd(j)  fscanf(entrada,"%D",&j)
#define EFSCAN(j)  fscanf(entrada,"%lf",&j)
#define EFSCANs(j) fscanf(entrada,"%s",j)

FILE *entrada;
char cod[10];

void load_g(char *Arq)
{
  extern double **g,xmin,xmax,ymin,ymax,gmin,gmax;
  extern long td_x,td_y;
  long cont1,cont2;

 /* Rotina LOADG para entrada de dados em formato .grd

    Parametros de Entrada:
      Nome do arquivo contendo os dados: arq

    Parametros de Saida:
      Ponteiro para o vetor das coordenadas em x da malha: x
      Idem para y
      Ponteiro de Ponteiros para a matriz de dados: g
      Quantidade de dados em x e y: td_x e td_y

  */

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

 /* Alocacao de memoria */

 g = dmatrix((long)1,td_x,(long)1,td_y);

 for(cont2=1;cont2<(td_y+1);cont2++)
  for(cont1=1;cont1<(td_x+1);cont1++)
   EFSCAN(g[cont1][cont2]);
   
 fclose(entrada);

}
