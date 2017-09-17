#include <stdio.h>
#include <math.h>
#include <time.h>
#include "inverft.h"

#define PI (float)3.141592653
#define G (float)6.67259e-3

#define EFSCANd(j)  fscanf(entrada,"%D",&j)
#define EFSCAN(j)  fscanf(entrada,"%f",&j)
#define CFSCAN(j)  fscanf(constantes,"%f",&j)
#define CFSCANd(j)  fscanf(constantes,"%D",&j)
#define EFSCANs(j) fscanf(entrada,"%s",j)
#define SFPRINs(j) fprintf(saida,"%s",j)
#define SFPRINg(j) fprintf(saida,"%.8lf\t",j)
#define SFPRINd(j) fprintf(saida,"%d\t",j)
#define SCAN scanf("%D",&w)

FILE *constantes;
FILE *entrada;
FILE *saida;

void main(void);

char arquivo0[20]="inverft.dat"; 
char arquivo1[20]="input.grd"; 
char arquivo2[20]="outputft.grd";
char arq3[20]="rmsft.dat";   
char cod[10];      

float rho,zzero;
float ***Dg,**specDg;
float ***h;
float xmin,xmax;
float ymin,ymax,hmin,hmax;
float interv_x,interv_y;
long tdat_x,tdat_y;
float ***hn,***Hn,**sp,***Soma,**sm;
long n,fal;
long cont,cont1,cont2;
float norm_Hn,norm_sp,norm_Soma,norm_sm,norm_h;
long iter;
long w;
float dgmin,dgmax;
long Conv_Exit,Conv_N;
long C_Exit,C_N;
float DU,DV,WH,SH,Umax,Vmax,Kmax;
float pcentW,pcentS;
float aux1,aux2,*rmsh;
float razao,S2,Sn;
float delta;
long MAX_ITER;
time_t t1,t2;


/* INICIO DO PROGRAMA */

void main()
{
 t1 = time(NULL);
 
 if ((constantes = fopen(arquivo0,"r")) == NULL)
  {
   fprintf(stderr,"Erro na abertura do arquivo %s",arquivo0);
   exit(-1);
  }
		     
 CFSCAN(rho);         /* contraste de densidades [g/cm3] */

 CFSCAN(zzero);       /* profundidade do nivel zero p/ topografia [metros] */
 CFSCAN(pcentW);      /* inicio do filtro passa baixa  [% da maior frequencia] (12.) */
 CFSCAN(pcentS);      /* final do filtro [idem] (12.5) */
 CFSCAN(delta);       /* convergencia para a razao Sn/S2 [parametro de inversao (0.001) */
 CFSCANd(MAX_ITER);   /* numero maximo de iteracoes */

 rmsh = vector((long)2,MAX_ITER);
 
 if ((entrada = fopen(arquivo1,"r")) == NULL)
  {
   fprintf(stderr,"Erro na abertura do arquivo %s",arquivo1);
   exit(-1);
  }

 EFSCANs(cod);
 EFSCANd(tdat_x); EFSCANd(tdat_y);
 EFSCAN(xmin); EFSCAN(xmax);
 EFSCAN(ymin); EFSCAN(ymax);
 EFSCAN(dgmin); EFSCAN(dgmax);

 Dg = f3tensor((long)1,(long)1,(long)1,tdat_x,(long)1,tdat_y);
 specDg = matrix((long)1,(long)1,(long)1,2 * tdat_y);
 h = f3tensor((long)1,(long)1,(long)1,tdat_x,(long)1,tdat_y);

 for(cont2=1;cont2<(tdat_y+1);cont2++)
   for(cont1=1;cont1<(tdat_x+1);cont1++)
     EFSCAN(Dg[1][cont1][cont2]);

 /* Soma de um plato de Bouguer com a espessura = zzero */

 for(cont2=1;cont2<(tdat_y+1);cont2++)
  for(cont1=1;cont1<(tdat_x+1);cont1++)
    Dg[1][cont1][cont2] = Dg[1][cont1][cont2] - ( (double)2. * PI * G * rho * zzero);

 /* FFT da matriz 'Dg' e sua multiplicacao por funcao do no.de onda */

 for (cont1=1;cont1<(2 * tdat_y + 1);cont1++)
   specDg[1][cont1] = (float)0.;

 (void) rlft3(Dg,specDg,(unsigned long)1,(unsigned long)tdat_x,(unsigned long)tdat_y,(int)1);

 interv_x = (xmax - xmin)/(float)(tdat_x - 1);
 interv_y = (ymax - ymin)/(float)(tdat_y - 1);
 (void) Gerokzi(Dg,specDg,tdat_x,tdat_y,interv_x,interv_y);

 /* Calculo de elementos do filtro passa-baixas */

 DU = ((float)2. * PI / (((float)tdat_x)*interv_x));
 DV = ((float)2. * PI / (((float)tdat_y)*interv_y));
 Umax = (float)(tdat_x/(long)2) * DU;
 Vmax = (float)(tdat_y/(long)2) * DV;
 Kmax = (float)sqrt((double)(Umax*Umax + Vmax*Vmax));
 WH = Kmax * (pcentW / (float)100.);
 SH = Kmax * (pcentS / (float)100.);

 /* Aplicacao do filtro passa-baixas a anomalia gravimetrica */
 (void) Filtro(Dg,specDg,tdat_x,tdat_y,interv_x,interv_y,WH,SH);

 /* Alocacao de memoria para elementos de trabalho */

 hn = f3tensor((long)1,(long)1,(long)1,tdat_x,(long)1,tdat_y);
 Hn = f3tensor((long)1,(long)1,(long)1,tdat_x,(long)1,tdat_y);
 sp = matrix((long)1,(long)1,(long)1,2 * tdat_y);
 Soma = f3tensor((long)1,(long)1,(long)1,tdat_x,(long)1,tdat_y);
 sm = matrix((long)1,(long)1,(long)1,2 * tdat_y);

 /* Estimativa inicial da topografia 'h()': zeros */
 for (cont1=1;cont1<(tdat_x + 1);cont1++)
  for (cont2=1;cont2<(tdat_y + 1);cont2++)
    h[1][cont1][cont2] = (float)0.;
 
 /* INICIO DO LOOP DO PROCESSO ITERATIVO DE INVERSAO  */

 printf("Iteracao n.01: topografia de zeros");

 for(iter=2;iter<(MAX_ITER + 1);iter++)
 {
  printf("\n");
  printf("Iteracao n. %d\n",iter);

 
  /* Preenchimento da matriz hn com 'hs'
     Alimentacao das matrizes Soma e sm com
     zeros. */

  for (cont1=1;cont1<(tdat_x + 1);cont1++)
   for (cont2=1;cont2<(tdat_y + 1);cont2++)
     hn[1][cont1][cont2] = h[1][cont1][cont2];

  for (cont1=1;cont1<(tdat_x + 1);cont1++)
   for (cont2=1;cont2<(tdat_y + 1);cont2++)
     Soma[1][cont1][cont2] = (float)0.;

  for (cont1=1;cont1<(2 * tdat_y + 1);cont1++)
   sm[1][cont1] = (float)0.;

  /* INICIO DO LOOP DE INCREMENTO DA SOMA <- espectro de hn vezes parcel */
  Conv_Exit = (long)0;
  Conv_N = (long)1;

  while (Conv_Exit != (long)1.)
   {
    Conv_N++;
    /* Calculo de Fatorial */
    fal = (long)1;
    for(cont=1;cont<(Conv_N+1);cont++) fal = fal * cont;

    /* Alimentacao da matriz hn com o produto h*hn (termo a termo) */
    for (cont1=1;cont1<(tdat_x + 1);cont1++)
      for (cont2=1;cont2<(tdat_y + 1);cont2++)
	hn[1][cont1][cont2] = hn[1][cont1][cont2] * h[1][cont1][cont2];

    /* Copia da matriz hn em Hn */
    for (cont1=1;cont1<(tdat_x + 1);cont1++)
      for (cont2=1;cont2<(tdat_y + 1);cont2++)
	Hn[1][cont1][cont2] = hn[1][cont1][cont2];

    /* Alimentacao da matriz auxiliar 'sp' com zeros */
    for (cont1=1;cont1<(2 * tdat_y + 1);cont1++)
      sp[1][cont1] = (float)0.;

    /* Aplicacao da Transformada de Fourier em 2 dimensoes */
    (void) rlft3(Hn,sp,(unsigned long)1,(unsigned long)tdat_x,(unsigned long)tdat_y,(int)1);

    /* Multiplicacao de Hn e sp por uma funcao do n.de onda */
    (void) parcel(Hn,sp,tdat_x,tdat_y,interv_x,interv_y,fal,Conv_N);

    /* Aplicacao do filtro passa-baixas ao espectro da topografia */
    (void) Filtro(Hn,sp,tdat_x,tdat_y,interv_x,interv_y,WH,SH);

    /* APLICACAO DO CRITERIO DE CONVERGENCIA - Espectro da topografia */

    if (iter == (long)2) Conv_Exit=(long)1;
    else
     {
      if (Conv_N == (long)2)
       {
	/* Multiplicacao de Hn e sp por uma funcao do n.de onda
	   e calculo do maximo do modulo dos elem. da matriz   */

	printf("Calculo de S2\n");
	S2 = maxover(Hn,sp,tdat_x,tdat_y,interv_x,interv_y);
	printf("S2 = %.8f\n",S2);
	razao = (float)1.;
       }
      else
       {
	/* Multiplicacao de Hn e sp por uma funcao do n.de onda
	 e calculo do maximo do modulo dos elem. da matriz     */

	Sn = maxover(Hn,sp,tdat_x,tdat_y,interv_x,interv_y);
	razao = Sn/S2;
	printf("razao = %.8f\n",razao);
       }

      if (razao < delta) Conv_Exit = (long)1;
      else
       {
	(void) Perfaz_Soma1(Soma,Hn,tdat_x,tdat_y);
	(void) Perfaz_Soma2(sm, sp, tdat_y);
       }
     }

   } /* fim do loop em Conv_Exit */

  printf("Espectro: convergiu em %d iteracoes\n",(Conv_N-(long)1));

  /* Subtracao das matrizes Dg e specDg pelas matrizes Soma e sm */

  (void) Perfaz_Subt1(Soma,Dg,tdat_x,tdat_y);
  (void) Perfaz_Subt2(sm,specDg, tdat_y);

  /* Aplicacao da transformada inversa para obtencao da proxima
     estimativa para a topografia.                             */

  (void) rlft3(Soma,sm,(unsigned long)1,(unsigned long)tdat_x,(unsigned long)tdat_y,(int)(-1));

  /* Multiplicacao de uma constante (2/(1*tdat_y*tdat_y)) */

  aux1 =  ((float)2.) / ( ((float)tdat_x) * ((float)tdat_y) );
  for (cont1=1;cont1<(tdat_x + 1);cont1++)
    for (cont2=1;cont2<(tdat_y + 1);cont2++)
       Soma[1][cont1][cont2] = aux1 * Soma[1][cont1][cont2];

  /* APLICACAO DE CRITERIO DE CONVERGENCIA - Estimativa de h */

  if (iter == (long) 2)
   {
    aux1 = (float) 0.;
    for (cont1=1;cont1<(tdat_x + 1);cont1++)
      for (cont2=1;cont2<(tdat_y + 1);cont2++)
	aux1 = aux1 + (Soma[1][cont1][cont2] * Soma[1][cont1][cont2]);

    rmsh[iter] = (float)sqrt((double) (aux1 /(tdat_x * tdat_y)));
   }
  
  else
   {
    aux1 = (float) 0.;
    aux2 = (float) 0.;
    for (cont1=1; cont1 < (tdat_x + 1); cont1++)
      for (cont2=1; cont2 < (tdat_y + 1); cont2++)
       {
	aux1 = Soma[1][cont1][cont2] - h[1][cont1][cont2];
	aux2 = aux2 + (aux1 * aux1);
       }
    rmsh[iter] = (float) sqrt((double) (aux2  / (tdat_x * tdat_y)));
   }
  
  printf("RMS da diferenca hn+1 - hn : %.5f\n",rmsh[iter]);

  /*  NOVA ESTIMATIVA DA TOPOGRAFIA: h <- Soma  */

  for (cont1=1;cont1<(tdat_x + 1);cont1++)
    for (cont2=1;cont2<(tdat_y + 1);cont2++)
      h[1][cont1][cont2] = Soma[1][cont1][cont2];
 
 } /* FINAL DO LOOP PRINCIPAL */

 (void) saverms(arq3, rmsh);

 
 /* Salvando a topografia em arquivo */

 (void) mudaorizf(h,zzero);

 hmin = h[1][1][1];
 hmax = h[1][1][1];

 for (cont1=1;cont1<(tdat_x + 1);cont1++)
   for (cont2=1;cont2<(tdat_y + 1);cont2++)
     {
      if (hmin > h[1][cont1][cont2])
	hmin = h[1][cont1][cont2];
      if (hmax < h[1][cont1][cont2])
	hmax = h[1][cont1][cont2];
     }

 if ((saida = fopen(arquivo2,"w")) == NULL)
  {
   fprintf(stderr,"Erro na abertura do arquivo %s",arquivo2);
   exit(-1);
  }

 SFPRINs("DSAA\n");
 SFPRINd(tdat_x); SFPRINd(tdat_y);
 SFPRINs("\n");
 SFPRINg((float)xmin); SFPRINg((float)xmax);
 SFPRINs("\n");
 SFPRINg((float)ymin); SFPRINg((float)ymax);
 SFPRINs("\n");
 SFPRINg((float)hmin); SFPRINg((float)hmax);
 SFPRINs("\n");

 for(cont2=1;cont2<(tdat_y+1);cont2++)
  {
   for(cont1=1;cont1<(tdat_x+1);cont1++) SFPRINg(h[1][cont1][cont2]);
   SFPRINs("\n");
  }

 /* Livrando espaco na memoria */

 (void) free_f3tensor(h,(long)1,(long)1,(long)1,tdat_x,(long)1,tdat_y);
 (void) free_f3tensor(hn,(long)1,(long)1,(long)1,tdat_x,(long)1,tdat_y);
 (void) free_f3tensor(Dg,(long)1,(long)1,(long)1,tdat_x,(long)1,tdat_y);
 (void) free_f3tensor(Hn,(long)1,(long)1,(long)1,tdat_x,(long)1,tdat_y);
 (void) free_matrix(sp,(long)1,(long)1,(long)1,2 * tdat_y);
 (void) free_matrix(specDg,(long)1,(long)1,(long)1,2 * tdat_y);
 (void) free_f3tensor(Soma,(long)1,(long)1,(long)1,tdat_x,(long)1,tdat_y);
 (void) free_matrix(sm,(long)1,(long)1,(long)1,2 * tdat_y);

 t2 = time(NULL);

 printf("Tempo decorrido [segundos]: %d\n",t2-t1);


} /* Fim do Programa */




void parcel(float ***H,float **spec,long N,long M,float Dx,float Dy,long f,long n0)
 {
  float u,v;
  float Dv,Du,fator;
  long i,j;

  Du = ((float)2. * PI / (((float)N)*Dx));
  Dv = ((float)2. * PI / (((float)M)*Dy));

  for(j=1;j<((M/2) +1);j++)
   {
    v = (float)(j - 1) * Dv;
    for(i=1;i<((N/2) +1);i++)
     {
      u = (float)(i - 1) * Du;
      fator = ffator(u,v,f,n0);
      H[1][i][2*j -1] = fator*H[1][i][2*j -1];
      H[1][i][2*j]    = fator*H[1][i][2*j];
     }
    for(i=((N/2) +1);i<(N +1);i++)
     {
      u = (float)(i - N - 1) * Du;
      fator = ffator(u,v,f,n0);
      H[1][i][2*j -1] = fator*H[1][i][2*j -1];
      H[1][i][2*j]    = fator*H[1][i][2*j];
     }
   }

  v = (float)(M/2) * Dv;
    for(i=1;i<((N/2) +1);i++)
     {
      u = (float)(i - 1) * Du;
      fator = ffator(u,v,f,n0);
      spec[1][2 * i - 1] = fator*spec[1][2 * i - 1];
      spec[1][2 * i] = fator*spec[1][2 * i];
     }
    for(i=((N/2) +1);i<(N +1);i++)
     {
      u = (float)(i - N - 1) * Du;
      fator = ffator(u,v,f,n0);
      spec[1][2 * i - 1] = fator*spec[1][2 * i - 1];
      spec[1][2 * i] = fator*spec[1][2 * i];
     }
 }


float ffator(float U, float V, long ffal, long nn)
 {
  float k,potk,resul;

  k = (float)sqrt((double)(U*U + V*V));
  if (nn==(long)1)
    potk = (float)1.;
  else
    potk = (float)pow((double)k,(double)(nn-1));
  resul = potk/(float)ffal;

  return(resul);
 }


void Perfaz_Soma1(float ***A, float ***B, long N, long M)
 {
  long ca,cb;

  for (ca=1;ca<(N + (long)1);ca++)
    for (cb=1;cb<(M + (long)1);cb++)
      A[1][ca][cb] = A[1][ca][cb] + B[1][ca][cb];
 }

void Perfaz_Soma2(float **a, float **b, long M)
 {
  long ca;
  for (ca=1;ca<(2 * M + (long)1);ca++)
    a[1][ca] = a[1][ca] + b[1][ca];
 }

void Perfaz_Subt1(float ***A, float ***B, long N, long M)
 {
  long ca,cb;

  for (ca=1;ca<(N + (long)1);ca++)
    for (cb=1;cb<(M + (long)1);cb++)
      A[1][ca][cb] =  (float)(-1.) * ( B[1][ca][cb] + A[1][ca][cb] );
 }

void Perfaz_Subt2(float **a, float **b, long M)
 {
  long ca;
  for (ca=1;ca<(2 * M + (long)1);ca++)
    a[1][ca] =  (float)(-1.) * ( b[1][ca] + a[1][ca] );
 }


void Gerokzi(float ***So,float **so,long N,long M,float Dx,float Dy)
 {
  float u,v;
  float Dv,Du,fator;
  long i,j;

  Du = ((float)2. * PI / (((float)N)*Dx));
  Dv = ((float)2. * PI / (((float)M)*Dy));

  for(j=1;j<((M/2) +1);j++)
   {
    v = (float)(j - 1) * Dv;
    for(i=1;i<((N/2) +1);i++)
     {
      u = (float)(i - 1) * Du;
      fator = Fator2(u,v,N,M);
      So[1][i][2*j -1] = fator*So[1][i][2*j -1];
      So[1][i][2*j]    = fator*So[1][i][2*j];
     }
    for(i=((N/2) +1);i<(N +1);i++)
     {
      u = (float)(i - N - 1) * Du;
      fator = Fator2(u,v,N,M);
      So[1][i][2*j -1] = fator*So[1][i][2*j -1];
      So[1][i][2*j]    = fator*So[1][i][2*j];
     }
   }

  v = (float)(M/2) * Dv;
    for(i=1;i<((N/2) +1);i++)
     {
      u = (float)(i - 1) * Du;
      fator = Fator2(u,v,N,M);
      so[1][2 * i - 1] = fator*so[1][2 * i - 1];
      so[1][2 * i] = fator*so[1][2 * i];
     }
    for(i=((N/2) +1);i<(N +1);i++)
     {
      u = (float)(i - N - 1) * Du;
      fator = Fator2(u,v,N,M);
      so[1][2 * i - 1] = fator*so[1][2 * i - 1];
      so[1][2 * i] = fator*so[1][2 * i];
     }
 }


float Fator2(float U, float V, long N, long M)
 {
  float k,resul,esc;
  k = (float)sqrt((double)(U*U + V*V));
  esc = (float)1.;
  resul = ( esc / ( (float)(2.) * PI * G * rho) ) * (float)exp( (double)(k * zzero) );
  return(resul);
 }


float maxover(float ***So,float **so,long N,long M,float Dx,float Dy)
 {
  float u,v;
  float Dv,Du,fator,aux;
  long i,j;
  float cand_max,Max=(float)0.;

  Du = ((float)2. * PI / (((float)N)*Dx));
  Dv = ((float)2. * PI / (((float)M)*Dy));

  for(j=1;j<((M/2) +1);j++)
   {
    v = (float)(j - 1) * Dv;
    for(i=1;i<((N/2) +1);i++)
     {
      u = (float)(i - 1) * Du;
      fator = FatorC(u,v);
      aux = fator * fator;
      cand_max = aux * (So[1][i][2*j -1] * So[1][i][2*j -1] + So[1][i][2*j] * So[1][i][2*j]);
      if (cand_max > Max) Max = cand_max;
     }
    for(i=((N/2) +1);i<(N +1);i++)
     {
      u = (float)(i - N - 1) * Du;
      fator = FatorC(u,v);
      aux = fator * fator;
      cand_max = aux * (So[1][i][2*j -1] * So[1][i][2*j -1] + So[1][i][2*j] * So[1][i][2*j]);
      if (cand_max > Max) Max = cand_max;
     }
   }

  v = (float)(M/2) * Dv;
    for(i=1;i<((N/2) +1);i++)
     {
      u = (float)(i - 1) * Du;
      fator = FatorC(u,v);
      aux = fator * fator;
      cand_max = aux * (so[1][2*i -1] * so[1][2*i -1] + so[1][2*i] * so[1][2*i]);
      if (cand_max > Max) Max = cand_max;
     }
    for(i=((N/2) +1);i<(N +1);i++)
     {
      u = (float)(i - N - 1) * Du;
      fator = FatorC(u,v);
      aux = fator * fator;
      cand_max = aux * (so[1][2*i -1] * so[1][2*i -1] + so[1][2*i] * so[1][2*i]);
      if (cand_max > Max) Max = cand_max;
     }

   Max = (float)sqrt((double)Max);
   return(Max);
 }

float FatorC(float U, float V)
 {
  float k,resul,esc;
  k = (float)sqrt((double)(U*U + V*V));
  resul = (float)exp((double)((float)(-1.) * k * zzero));
  return(resul);
 }


void Filtro(float ***H,float **spec,long N,long M,float Dx,float Dy,float wh, float sh)
 {
  float u,v;
  float Dv,Du,fator;
  long i,j;

  Du = ((float)2. * PI / (((float)N)*Dx));
  Dv = ((float)2. * PI / (((float)M)*Dy));

  for(j=1;j<((M/2) +1);j++)
   {
    v = (float)(j - 1) * Dv;
    for(i=1;i<((N/2) +1);i++)
     {
      u = (float)(i - 1) * Du;
      fator = filtro(u,v,wh,sh);
      H[1][i][2*j -1] = fator*H[1][i][2*j -1];
      H[1][i][2*j]    = fator*H[1][i][2*j];
     }
    for(i=((N/2) +1);i<(N +1);i++)
     {
      u = (float)(i - N - 1) * Du;
      fator = filtro(u,v,wh,sh);
      H[1][i][2*j -1] = fator*H[1][i][2*j -1];
      H[1][i][2*j]    = fator*H[1][i][2*j];
     }
   }

  v = (float)(M/2) * Dv;
    for(i=1;i<((N/2) +1);i++)
     {
      u = (float)(i - 1) * Du;
      fator = filtro(u,v,wh,sh);
      spec[1][2 * i - 1] = fator*spec[1][2 * i - 1];
      spec[1][2 * i] = fator*spec[1][2 * i];
     }
    for(i=((N/2) +1);i<(N +1);i++)
     {
      u = (float)(i - N - 1) * Du;
      fator = filtro(u,v,wh,sh);
      spec[1][2 * i - 1] = fator*spec[1][2 * i - 1];
      spec[1][2 * i] = fator*spec[1][2 * i];
     }
 }


float filtro(float U, float V, float ww, float ss)
 {
  float k,resul,aux1,aux2;
  k = (float)sqrt((double)(U*U + V*V));
  if (k<ww)
    resul = (float)1.;
  else
    if (k>ss)
      resul = (float)0.;
    else
      {
       aux1 = PI * ((k - ww)/(ss - ww));
       aux2 = (float)cos( (double)aux1 );
       resul = (float)0.5 * ( (float)1. + aux2 );
      }
  return(resul);
 }


