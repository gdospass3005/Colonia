#include<stdio.h>
#include<math.h>
#include"prisboug.h"

#define PRINd(j)  fprintf(saidarms,"%d\t",j)
#define PRIN(j)  fprintf(saidarms,"%.8G\t",j)
#define PRINs(j)  fprintf(saidarms,"%s",j)

FILE *saidarms;

void saverms(char *Arq3, double *RMSH, double *RMSG, double *MAXERRG, double *MAXERRH)
{

  extern long MAX_ITER;

  long cont1;

 if ((saidarms = fopen(Arq3,"w")) == NULL)
  {
   fprintf(stderr,"Erro na abertura do arquivo %s",Arq3);
   exit(-1);
  }
 
 for(cont1=1;cont1<(MAX_ITER+1);cont1++)
  {
   PRINd(cont1);
   PRIN(RMSH[cont1]); PRIN(MAXERRH[cont1]); PRIN(RMSG[cont1]); PRIN(MAXERRG[cont1]);
   PRINs("\n");
  }

 fclose(saidarms);
}
