#include<stdio.h>
#include<math.h>
#include"prisboug.h"

#define PI (float)3.141592653

void loop(double **G, double **H)
{
 extern double *rmsh,*rmsg,*maxerrg,*maxerrh,xmin,xmax,ymin,ymax,HMAX;
 extern long td_x,td_y,MAX_ITER;
 extern double a0,a1,a2;

 double i_antx,i_posx,i_anty,i_posy,peso;
 double h0_antx,h0_posx,h0_anty,h0_posy,h0;
 double gamma=(double)6.67259e-3; 
 double R2,dz;
 long cont1,cont2,iter;
 double **gg,i,j,k,m,n,p,mg,mh;

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
	h0_antx = H[cont1-1][cont2];
	h0_posx = H[cont1+1][cont2];
	h0_anty = H[cont1][cont2-1];
	h0_posy = H[cont1][cont2+1];
	
	R2 = a0 + a1 * h0 + a2 * h0 * h0;
	dz=(G[cont1][cont2]-gg[cont1][cont2])/((double)2. * PI * gamma * R2);
	i = h0 + dz;

	if((cont1>1) && (cont1<td_x) && (cont2>1) && (cont2<td_y))
	 {
	  R2=a0 +a1 * h0_antx +a2 * h0_antx * h0_antx;
	  dz=(G[cont1-1][cont2]-gg[cont1-1][cont2])/((double)2. * PI * gamma * R2);
	  i_antx = h0_antx + dz;
	  
	  R2=a0 +a1 * h0_posx +a2 * h0_posx * h0_posx;
	  dz=(G[cont1+1][cont2]-gg[cont1+1][cont2])/((double)2. * PI * gamma * R2);
	  i_posx = h0_posx + dz;
	
	  R2=a0 +a1 * h0_anty +a2 * h0_anty * h0_anty;
	  dz=(G[cont1][cont2-1]-gg[cont1][cont2-1])/((double)2. * PI * gamma * R2);
	  i_anty = h0_anty + dz;

	  R2=a0 +a1 * h0_posy +a2 * h0_posy * h0_posy;
	  dz=(G[cont1][cont2+1]-gg[cont1][cont2+1])/((double)2. * PI * gamma * R2);
	  i_posy = h0_posy + dz;
	  
	  peso = fabs(h0 / HMAX);
	  
	  i = (double)0.25 * peso * (i_antx + i_posx + i_anty + i_posy) + ((double)1. - peso) * i;
	  
	 } 
	 
	if(i > HMAX) i = HMAX;
	if(i < ((double)-1. * HMAX)) i = -HMAX;

	m = i - h0;     /* Variacao da topografia (n,n-1) */  
	j = j + (m * m);             /* Somatoria quadratica */

	n = fabs(G[cont1][cont2]-gg[cont1][cont2]);  /* Difs. de gobs - gcal */
 
	if(n >= mg) maxerrg[iter-1] = n; /* Maximo residuo da anomalia (Cordell) */
	else maxerrg[iter-1] = mg;
 
	mg = maxerrg[iter-1];

	p = fabs(h0-i);  /* Variacao de Hn - Hn-1 */
 
	if(p >= mh) maxerrh[iter] = p; /* Maxima variacao da topografia */
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
