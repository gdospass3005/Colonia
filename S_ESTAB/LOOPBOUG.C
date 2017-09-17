#include<stdio.h>
#include<math.h>
#include"prisver.h"

#define PI (float)3.141592653

void loop(double **G, double **H)
{
 extern double *rmsh,*rmsg,*maxerrg,*maxerrh,xmin,xmax,ymin,ymax,HMAX;
 extern long td_x,td_y,MAX_ITER;
 extern double a0,a1,a2;

 double h0;
 double gamma=(double)6.67259e-3; 
 double R2,dz;
 long cont1,cont2,iter;
 double **gg,i,j,k,m,p,n,mg,mh;

 gg = dmatrix((long)1,td_x,(long)1,td_y); 

 for(iter=1;iter<(MAX_ITER + 1);iter++)   /* As iteracoes comecam aqui */
  {
     j = (double)0.;
     k = (double)0.;
     mg = (double)0.; 
     mh = (double)0.; 

    (void) avalia(H, gg, iter);

    for(cont2=1;cont2<(td_y+1);cont2++)
     for(cont1=1;cont1<(td_x+1);cont1++)
       {
	h0 = H[cont1][cont2];
	
	R2 = a0 + a1 * h0 + a2 * h0 * h0;
	dz=(G[cont1][cont2]-gg[cont1][cont2])/((double)2. * PI * gamma * R2);
	i = h0 + dz;

	m = i - h0;     /* Variacao da topografia (n,n-1) */  
	j = j + (m * m);             /* Somatoria quadratica */

	n = fabs(G[cont1][cont2]-gg[cont1][cont2]);  /* Difs. de gobs - gcal */
 
	if(n >= mg) maxerrg[iter-1] = n; /* Maximo residuo da anomalia (Cordell) */
	else maxerrg[iter-1] = mg;
 
	mg = maxerrg[iter-1];

	p=fabs(h0-i);  /* Variacao Hn - Hn-1 */
	if(p >= mh) maxerrh[iter] = p; /* Maxima variacao Hn - Hn-1 */
	else maxerrh[iter] = mh;

	mh = maxerrh[iter];
	
	k = k + n * n;             /* Somatoria quadratica */
	H[cont1][cont2] = i;
       }
    
    rmsh[iter] = sqrt( j / ( td_x * td_y) );
    rmsg[iter-1] = sqrt( k / ( td_x * td_y) );
    
    printf(".");

  }
 printf("\n");

 (void) avalia(H, gg, (MAX_ITER+1)); /* Anomalia gerada pelo resultado H */
 k=(double)0.;
 mg=(double)0.;
 for(cont2=1;cont2<(td_y+1);cont2++)
  for(cont1=1;cont1<(td_x+1);cont1++)
    {
     n=fabs(G[cont1][cont2]-gg[cont1][cont2]);  /* Difs. de gobs - gcal */
       
     if(n >= mg) maxerrg[MAX_ITER] = n; /* Maximo residuo da anomalia (Cordell) */
     else maxerrg[MAX_ITER] = mg;
     mg = maxerrg[MAX_ITER];
     k = k + n * n;         /* Somatoria quadratica */
     G[cont1][cont2] = gg[cont1][cont2];
    }
 rmsg[MAX_ITER] = sqrt( k / ( td_x * td_y) );

 (void) free_dmatrix(gg,(long)1,td_x,(long)1,td_y);
} 
