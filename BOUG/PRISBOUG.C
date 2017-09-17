/* Programa PRISBOUG

   Avalia profundidade das bases de prismas 3D a partir de anomalia grav.
   (mais fiel ao programa de Rao & Babu, 1990)

   Entrada: Grid com Anomalia - input.grd
	    Coeficientes da funcao quadratica de densidade - prisver.dat
	    Valores de MAX_ITER,ITER1,ITER2,HMAX

   Saida:  Grid com profundidades da base dos prismas      - outputpv.grd
	   Dados de RMS:                                   - rmspv.dat
	   N.da iteracao, rms das topografias(*), rms das anomalias(**)
	   (*): anterior menos atual
	   (**): observada menos calculada

*/         

#include <math.h>
#include <stdio.h>
#include <time.h>
#include "prisboug.h"

void main(void);

char arq[20] =  "input.grd";
char arq1[20] = "prisboug.dat";
char arq2[20] = "outputbo.grd";
char arq3[20] = "rmsboug.dat";
char arq4[20] = "anomalbo.grd";

double **h,**g,*rmsh,*rmsg,*maxerrg,*maxerrh,xmin,xmax,ymin,ymax,gmin,gmax;
double interv_x,interv_y,difinterv,dx,dy,t,w;
long td_x,td_y;
long LT1=(long)1;
long LT2=(long)3;
long LT3=(long)5;

double a0,a1,a2;    /* Coeficiente da funcao contraste de densidades  */

long MAX_ITER;      /* Limite maximo para a quantidade de iteracoes   */

long ITER1,ITER2;   /* ITER1 = 1 vizinho, ITER2 = 3 vizinhos,         */
		    /* Da iteracao (ITER2 + 1) em diante = 5 vizinhos */
		    /* (para a utilizacao da eq.exata p/ a anomalia   */
		    /*  prisma dado).                                 */

double HMAX; /* Maxima profundidade estimada para as fontes 
		(base dos prismas) <--> parametro de estabilizacao */
	     
time_t t1,t2;

void main()
{
 t1 = time(NULL);

 (void) load_g(arq);
 (void) load_coef(arq1);

 h = dmatrix((long)1,td_x,(long)1,td_y);
 rmsh = dvector((long)1, MAX_ITER ); 
 rmsg = dvector((long)0, MAX_ITER ); 
 maxerrg = dvector((long)0, MAX_ITER );
 maxerrh = dvector((long)1, MAX_ITER );

 (void) flat(g, h);
 (void) loop(g, h); 
 (void) save(arq2, h);
 (void) saverms(arq3, rmsh, rmsg, maxerrg, maxerrh);
 (void) save(arq4, g);


 (void) free_dmatrix(g,(long)1,td_x,(long)1,td_y);
 (void) free_dmatrix(h,(long)1,td_x,(long)1,td_y);
 (void) free_dvector(rmsh,(long)1, MAX_ITER);
 (void) free_dvector(rmsg,(long)0, MAX_ITER);
 (void) free_dvector(maxerrg,(long)0, MAX_ITER);
 (void) free_dvector(maxerrh,(long)1, MAX_ITER);

 t2 = time(NULL);
 t2 = t2 - t1;
 printf("Tempo decorrido [segundos]: %ld\n",t2);
}
