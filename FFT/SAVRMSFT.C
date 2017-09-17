#include<stdio.h>
#include<math.h>
#include"inverft.h"

#define PRINd(j)  fprintf(saidarms,"%d\t",j)
#define PRIN(j)  fprintf(saidarms,"%e\t",j)
#define PRINs(j)  fprintf(saidarms,"%s",j)

FILE *saidarms;

void saverms(char *Arq3, float *RMSH)
{ extern long MAX_ITER;
  long cont1;

 if ((saidarms = fopen(Arq3,"w")) == NULL)
  {
   fprintf(stderr,"Erro na abertura do arquivo %s",Arq3);
   exit(-1);
  }
 
 for(cont1=2;cont1<(MAX_ITER+1);cont1++)
  {
   PRINd(cont1);
   PRIN(RMSH[cont1]);
   PRINs("\n");
  }

 fclose(saidarms);
}
