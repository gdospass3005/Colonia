#include<stdio.h>
#include<math.h>
#include"prisft.h"

void loopft(double **G, double **H)
{

 extern double *rmsh,*rmsg,*maxerrg,*maxerrh,xmin,xmax,ymin,ymax,HMAX;
 extern long td_x,td_y,MAX_ITER;

 double i_antx,i_posx,i_anty,i_posy,peso,sgn,mg,mh;
 
 long cont1,cont2,iter;
 double **gg,i,j,k,m,n,p;

 gg = dmatrix((long)1,td_x,(long)1,td_y); 

 for(iter=1;iter<(MAX_ITER + 1);iter++)   /* As iteracoes comecam aqui */
  {
     j = (double)0.;
     k = (double)0.;
     mg=(double)0.;
     mh=(double)0.;

    (void) avaliaft(H, gg);

    for(cont2=1;cont2<(td_y+1);cont2++)
     for(cont1=1;cont1<(td_x+1);cont1++)
       {
	sgn = H[cont1][cont2] / (fabs(H[cont1][cont2]));

	i = (G[cont1][cont2]/gg[cont1][cont2]) * H[cont1][cont2];

	if((cont1>1) && (cont1<td_x) && (cont2>1) && (cont2<td_y))
	 {
	  i_antx = (G[cont1 -1][cont2]/gg[cont1 -1][cont2]) * H[cont1 -1][cont2];
	  i_posx = (G[cont1 +1][cont2]/gg[cont1 +1][cont2]) * H[cont1 +1][cont2];
	  i_anty = (G[cont1][cont2 -1]/gg[cont1][cont2 -1]) * H[cont1][cont2 -1];
	  i_posy = (G[cont1][cont2 +1]/gg[cont1][cont2 +1]) * H[cont1][cont2 +1];

	  peso = sgn * (H[cont1][cont2] / HMAX);
	  
	  i = (double)0.25 * peso * (i_antx + i_posx + i_anty + i_posy) + ((double)1. - peso) * i;
	  
	 } 
	 
	if(i > HMAX) i = HMAX;
	if(i < ((double)-1. * HMAX)) i = -HMAX;

	m = i - H[cont1][cont2];    /* Difenrencas da topografia (n,n-1) */  
	j = j + (m * m);             /* Somatoria quadratica */
	n=fabs(G[cont1][cont2]-gg[cont1][cont2]);  /* Difs. de gobs - gcal */
	
	if(n >= mg) maxerrg[iter-1] = n; /* Maximo residuo da anomalia (Cordell) */
	else maxerrg[iter-1] = mg;
 
	mg = maxerrg[iter-1];

	p=fabs(H[cont1][cont2]-i);  /* Variacao Hn - Hn-1 */
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

 (void) avaliaft(H, gg); /* Anomalia gerada pelo resultado H */
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
