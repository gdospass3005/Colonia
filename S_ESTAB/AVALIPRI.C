#include<stdio.h>
#include<math.h>
#include <time.h>
#include "prisver.h"

static long dmaxarg1,dmaxarg2;
#define LMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ?\
	(dmaxarg1) : (dmaxarg2))

void avalia(double **H, double **G, long iter)
{
 extern double xmin,xmax,ymin,ymax,interv_x,interv_y,a0,a1,a2;
 extern long td_x,td_y;
 extern long LT1,LT2,LT3,ITER1,ITER2;

 long LT;
 long i,j,k,l,m,n,max_td,max_x,max_y,len_p,q;
 long l1,l2,l3,l4,delta_k,delta_l,limit2;
 double x0,y0,xp,yp;
 double z1=(double)0.000000001;
 double **gt;
 long cont1,cont2;
 double sgn;

 interv_x = (xmax - xmin)/(double)(td_x - 1);
 interv_y = (ymax - ymin)/(double)(td_y - 1);
 
 max_td = LMAX(td_x,td_y);
 gt = dmatrix((long)1,max_td,(long)1,max_td); /* Matriz para trapezios */

 if(iter <= ITER1) LT = LT1;
 if((iter > ITER1) && (iter <= (ITER2+ITER1))) LT = LT2;
 if(iter > (ITER2+ITER1)) LT = LT3;

 for(i=1;i<(td_x+1);i++)
  for(j=1;j<(td_y+1);j++)
   if( (H[i][j] < (double)0.00000001) && (H[i][j] > (double)(-0.00000001)) ) H[i][j] = (double)0.00000001;

 for(i=1;i<(td_x+1);i++)
   for(j=1;j<(td_y+1);j++)
     G[i][j]=(double)0.;

 x0 = (double)0.;
 y0 = x0;

 for(k=1;k<(td_x+1);k++)
  {
   for(l=1;l<(td_y+1);l++)
    {
     /* Calculo do sinal do prisma */
     sgn = H[k][l] / (fabs(H[k][l]));
     
     /* Calculo das dimensoes do trapezio de anomalias */

     l1 = k;
     l2 = td_x - k + (long)1;
     max_x = LMAX(l1,l2);
     l3 = l;
     l4 = td_y - l + (long)1;
     max_y = LMAX(l3,l4);
     len_p = LMAX(max_x,max_y);
     if (len_p == max_x) limit2 = max_y;
       else limit2 = max_x;

     /* Anomalias do trapezio ALTERNATIVA 'EXATA/APROX.'   */
     
     for(m=1;m<(limit2 + 1);m++)
      {
       if(len_p > (LT + 1))  
	{
	 for(n=m;n<(LT + 2);n++)
	  {
	   xp = ((double)(m-1))*interv_x;
	   yp = ((double)(n-1))*interv_y;
	   gt[m][n] =  anom1(x0,y0,xp,yp,z1,H[k][l]);
	  }
	 for(n=(LT + 2);n<(len_p + 1);n++)
	  {
	   xp = ((double)(m-1))*interv_x;
	   yp = ((double)(n-1))*interv_y;
	   gt[m][n] =  anom2(x0,y0,xp,yp,z1,H[k][l]);
	  }
	}
       else  
	 for(n=m;n<(len_p + 1);n++)
	  {
	   xp = ((double)(m-1))*interv_x;
	   yp = ((double)(n-1))*interv_y;
	   gt[m][n] =  anom1(x0,y0,xp,yp,z1,H[k][l]);
	  }
      }

     /* Contribuicao do prisma em cada ponto da malha */

     for(i=1;i<(td_x+1);i++)
      {
      for(j=1;j<(td_y+1);j++)
       {
	delta_k = (labs(k - i)) + (long)1;
	delta_l = (labs(l - j)) + (long)1;
	if ( delta_l >= delta_k )
	 {
	  m = delta_k;
	  n = delta_l;
	 }
	else
	 {
	  m = delta_l;
	  n = delta_k;
	 }
	G[i][j] = G[i][j] + (sgn * gt[m][n]);
       }
      }

    } /* Fim do loop em k */
  } /* Fim do loop em l */

(void) free_dmatrix(gt,(long)1,max_td,(long)1,max_td);

} /* Fim da rotina */
