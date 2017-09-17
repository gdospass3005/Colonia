/* Programa PRISFT --- Algoritmo Misto

   Avalia profundidade das bases de prismas 3D a partir de anomalia grav.
   gerada com a utilizacao da transformada de Fourier da interface
   (metodo de Parker).

   Entrada: Grid com Anomalia - input.grd
	    Coeficientes da funcao quadratica de densidade - inverft.dat

   Saida:  Grid com profundidades da base dos prismas - outputpf.grd
	   Dados de RMS                               - rmspf.dat
*/

#include <math.h>
#include <stdio.h>
#include <time.h>
#include "prisft.h"

void main(void);

char arq[20] =  "input.grd";
char arq1[20] = "prisft.dat";
char arq2[20] = "outputpf.grd";
char arq3[20] = "rmspf.dat";
char arq4[20] = "anomalpf.grd";

double **h,**g,*rmsh,*rmsg,*maxerrg,*maxerrh,xmin,xmax,ymin,ymax,gmin,gmax;
double interv_x,interv_y,difinterv;
long td_x,td_y;

double HMAX; /* Maxima profundidade estimada para as fontes 
		(base dos prismas) <--> parametro de estabilizacao */
	     

double rho,zzero,delta;
long MAX_ITER;

time_t t1,t2;

void main()
{

 t1 = time(NULL);

 (void) load_g(arq);

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

 (void) load_coef_ift(arq1);

 h = dmatrix((long)1,td_x,(long)1,td_y);
 rmsh = dvector((long)1, MAX_ITER ); 
 rmsg = dvector((long)0, MAX_ITER ); 
 maxerrg = dvector((long)0, MAX_ITER ); 
 maxerrh = dvector((long)1, MAX_ITER ); 

 (void) flat(g, h);

 (void) loopft(g, h); 

 (void) save(arq2, h);
 
 (void) saverms(arq3, rmsh, rmsg, maxerrg, maxerrh); 

 (void) save(arq4, g);

 (void) free_dmatrix(g,(long)1,td_x,(long)1,td_y);
 (void) free_dmatrix(h,(long)1,td_x,(long)1,td_y);

 (void) free_dvector(rmsg,(long) 0, MAX_ITER);
 (void) free_dvector(rmsh,(long) 1, MAX_ITER);
 (void) free_dvector(maxerrg,(long) 0, MAX_ITER);
 (void) free_dvector(maxerrh,(long) 1, MAX_ITER);

 t2 = time(NULL);
   t2 = t2 - t1;
   printf("Tempo decorrido [segundos]: %ld\n",t2);
}
