#include<stdio.h>
#include<math.h>
#include"prisboug.h"

#define PRINd(j)  fprintf(saida,"%d\t",j)
#define PRIN(j)  fprintf(saida,"%.8lf\t",j)
#define PRINs(j) fprintf(saida,"%s",j)

FILE *saida;

void save(char *Arq2, double **G)
{
  extern double xmin,xmax,ymin,ymax;
  extern long td_x,td_y;

  double gmax,gmin;
  long cont1,cont2;

/* Grava anomalia gravimetrica em arquivo formato ".grd"
   Rotina para o calculo de gmax, gmin. */

 gmin = G[1][1];
 gmax = G[1][1];

 for (cont1=1;cont1<(td_x + 1);cont1++)
   for (cont2=1;cont2<(td_y + 1);cont2++)
     {
      if (gmin > G[cont1][cont2])
	gmin = G[cont1][cont2];
      if (gmax < G[cont1][cont2])
	gmax = G[cont1][cont2];
     }

 if ((saida = fopen(Arq2,"w")) == NULL)
  {
   fprintf(stderr,"Erro na abertura do arquivo %s",Arq2);
   exit(-1);
  }

 PRINs("DSAA\n");
 PRINd(td_x); PRINd(td_y);
 PRINs("\n");
 PRIN(xmin); PRIN(xmax);
 PRINs("\n");
 PRIN(ymin); PRIN(ymax);
 PRINs("\n");
 PRIN(gmin); PRIN(gmax);
 PRINs("\n");

 for(cont2=1;cont2<(td_y+1);cont2++)
  {
   for(cont1=1;cont1<(td_x+1);cont1++) PRIN(G[cont1][cont2]);
   PRINs("\n");
  }

 fclose(saida);
}
