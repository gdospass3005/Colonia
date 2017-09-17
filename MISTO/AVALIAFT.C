#include<stdio.h>
#include<math.h>
#include "prisft.h"

#define PI (double)3.14159265358979
#define G (double)6.67259e-3  /* Unidades: g/cm3, metros e mGal  */

void avaliaft(double **H, double **gg)
{
 extern double interv_x,interv_y;
 extern double rho,zzero,delta;
 extern long td_x,td_y;

double **h,***hn,***Hn,**sp,***Soma,**sm;
long n,fal,tdat_x=td_x,tdat_y=td_y;
long cont,cont1,cont2;
double norm_Hn,norm_sp,norm_Soma,norm_sm,razao,S1,Sn;
double Bouguer;
long Conv_Exit,Conv_N;

/* Adaptacao para dados comuns:
   Mudancao na orientacao e origem das profundidades  */

   (void) mudaoriz(H,zzero);

/* Alocacao de memoria para elementos de trabalho */

 h = dmatrix((long)1,td_x,(long)1,td_y);
 
 for (cont1=1;cont1<(tdat_x + 1);cont1++)
   for (cont2=1;cont2<(tdat_y + 1);cont2++)
     h[cont1][cont2] = H[cont1][cont2];
 
 hn = d3tensor((long)1,(long)1,(long)1,td_x,(long)1,td_y);
 
 Hn = d3tensor((long)1,(long)1,(long)1,td_x,(long)1,td_y);
 
 sp = dmatrix((long)1,(long)1,(long)1,2 * td_y);
 Soma = d3tensor((long)1,(long)1,(long)1,td_x,(long)1,td_y);
 
 sm = dmatrix((long)1,(long)1,(long)1,2 * td_y);

/*  Preenchimento com '1s' da matriz hn
    Alimentacao das matrizes Soma, sm e sp com
    zeros. */

 for (cont1=1;cont1<(tdat_x + 1);cont1++)
   for (cont2=1;cont2<(tdat_y + 1);cont2++)
     hn[1][cont1][cont2] = (double)1.;

 for (cont1=1;cont1<(tdat_x + 1);cont1++)
   for (cont2=1;cont2<(tdat_y + 1);cont2++)
     Soma[1][cont1][cont2] = (double)0.;

 for (cont1=1;cont1<(2 * tdat_y + 1);cont1++)
   sm[1][cont1] = (double)0.;

/* INICIO DO LOOP DE INCREMENTO DA SOMA */
 Conv_Exit = (long)0.;
 Conv_N = (long)0.;

 while (Conv_Exit != (long)1.)
 {
 Conv_N++;
   /* printf("Iteracao n. :%d\n",Conv_N); */

/* Calculo de Fatorial */

 fal = (long)1;
 for(cont=1;cont<(Conv_N+1);cont++) fal = fal * cont;

 for (cont1=1;cont1<(2 * tdat_y + 1);cont1++)
   sp[1][cont1] = (double)0.;

/*  Alimentacao da matriz hn com o produto h*hn (termo a termo) */

 for (cont1=1;cont1<(tdat_x + 1);cont1++)
   for (cont2=1;cont2<(tdat_y + 1);cont2++)
     hn[1][cont1][cont2] = hn[1][cont1][cont2] * h[cont1][cont2];

/* Copia da matriz hn em Hn */

 for (cont1=1;cont1<(tdat_x + 1);cont1++)
   for (cont2=1;cont2<(tdat_y + 1);cont2++)
     Hn[1][cont1][cont2] = hn[1][cont1][cont2];

/* Aplicacao da Transformada de Fourier em 2 dimensoes */
 (void) rlft3(Hn,sp,(unsigned long)1,(unsigned long)tdat_x,(unsigned long)tdat_y,(int)1);

 /* printf(" Multiplicacao de Hn e sp por uma funcao do n.de onda \n"); */
 (void) parcel(Hn,sp,tdat_x,tdat_y,(double)interv_x,(double)interv_y,fal,Conv_N);


/* APLICACAO DO CRITERIO DE CONVERGENCIA */

 if (Conv_N==1)
  {
    /* Multiplicacao de Hn e sp por uma funcao do n.de onda
       e calculo do maximo do modulo dos elem. da matriz    */

    /* printf("Calculo de S1\n"); */
    S1 = maxover(Hn,sp,tdat_x,tdat_y,(double)interv_x,(double)interv_y);
    /* printf("S1 = %.8f\n",S1); */
    razao = (double)1.;
  }
 else
  {
    /* Multiplicacao de Hn e sp por uma funcao do n.de onda
       e calculo do maximo do modulo dos elem. da matriz    */

    Sn = (double)maxover(Hn,sp,tdat_x,tdat_y,(double)interv_x,(double)interv_y);
    /* printf("Sn = %.8f\n",Sn); */
    razao = Sn/S1;
    /* printf("razao = %.8f\n",razao); */
  }

 if (razao < delta) Conv_Exit = (long)1;
 else
  {
   (void) Perfaz_Soma1(Soma,Hn,tdat_x,tdat_y);
   (void) Perfaz_Soma2(sm, sp, tdat_y);
  }

 } /* fim do loop em Conv_Exit */

   /* printf("Convergiu em %d iteracoes\n",Conv_N); */

/* Etapa: multiplicacao da somatoria por um fator
   de escala que eh funcao de uma profundidade de
   referencia, do contraste de densidades e do
   n.de onda.                                    */

(void) Gerokz(Soma,sm,tdat_x,tdat_y,(double)interv_x,(double)interv_y);

/* Aplicacao da transformada inversa para obtencao da
   anomalia gravimetrica.                            */

 (void) rlft3(Soma,sm,(unsigned long)1,(unsigned long)tdat_x,(unsigned long)tdat_y,(int)(-1));

 for (cont1=1;cont1<(tdat_x + 1);cont1++)
   for (cont2=1;cont2<(tdat_y + 1);cont2++)
     gg[cont1][cont2] = (double)Soma[1][cont1][cont2];


 /* PLATO DE BOUGUER PARA A CORRECAO DE ZZERO    */

 Bouguer = (double)2. * PI * G * rho * zzero;
 for (cont1=1;cont1<(tdat_x + 1);cont1++)
   for (cont2=1;cont2<(tdat_y + 1);cont2++)
     gg[cont1][cont2] = gg[cont1][cont2] + Bouguer;


 (void) free_dmatrix(h,(long)1,(long)1,(long)1,2 * tdat_y);
 (void) free_d3tensor(hn,(long)1,(long)1,(long)1,tdat_x,(long)1,tdat_y);
 (void) free_d3tensor(Hn,(long)1,(long)1,(long)1,tdat_x,(long)1,tdat_y);
 (void) free_dmatrix(sp,(long)1,(long)1,(long)1,2 * tdat_y);
 (void) free_d3tensor(Soma,(long)1,(long)1,(long)1,tdat_x,(long)1,tdat_y);
 (void) free_dmatrix(sm,(long)1,(long)1,(long)1,2 * tdat_y);

  /*   Mudanca na orientacao e origem das profundidades  */

   (void) mudaoriz(H,zzero);

} /* Fim do Programa */


/* FUNCOES UTILIZADAS PELO PROGRAMA  */

void parcel(double ***H,double **spec,long N,long M,double Dx,double Dy,long f,long n0)
 {
  double u,v;
  double Dv,Du,fator;
  long i,j;

  Du = ((double)2. * PI / (((double)N)*Dx));
  Dv = ((double)2. * PI / (((double)M)*Dy));

  for(j=1;j<((M/2) +1);j++)
   {
    v = (double)(j - 1) * Dv;
    for(i=1;i<((N/2) +1);i++)
     {
      u = (double)(i - 1) * Du;
      fator = ffator(u,v,f,n0);
      H[1][i][2*j -1] = fator*H[1][i][2*j -1];
      H[1][i][2*j]    = fator*H[1][i][2*j];
     }
    for(i=((N/2) +1);i<(N +1);i++)
     {
      u = (double)(i - N - 1) * Du;
      fator = ffator(u,v,f,n0);
      H[1][i][2*j -1] = fator*H[1][i][2*j -1];
      H[1][i][2*j]    = fator*H[1][i][2*j];
     }
   }

  v = (double)(M/2) * Dv;
    for(i=1;i<((N/2) +1);i++)
     {
      u = (double)(i - 1) * Du;
      fator = ffator(u,v,f,n0);
      spec[1][2 * i - 1] = fator*spec[1][2 * i - 1];
      spec[1][2 * i] = fator*spec[1][2 * i];
     }

    for(i=((N/2) +1);i<(N +1);i++)
     {
      u = (double)(i - N - 1) * Du;
      fator = ffator(u,v,f,n0);
      spec[1][2 * i - 1] = fator*spec[1][2 * i - 1];
      spec[1][2 * i] = fator*spec[1][2 * i];
     }

 }

double ffator(double U, double V, long ffal, long nn)
 {
  double k,potk,resul;

  k = (double)sqrt((double)(U*U + V*V));
  if (nn==(long)1)
    potk = (double)1.;
  else
    potk = (double)pow((double)k,(double)(nn-1));

  resul = potk/ffal;

  return(resul);
 }

double norma_t(double ***A, long N, long M)
 {
  long c1,c2;
  double resul=(double)0.;

  for (c1=1;c1<(N + 1);c1++)
    for (c2=1;c2<(M + 1);c2++)
      resul = resul + A[1][c1][c2] * A[1][c1][c2];

  resul = (double)sqrt((double)resul);
  return(resul);
 }

double norma_m(double **A,long M)
 {
  long c1;
  double resul=(double)0.;

  for (c1=1;c1<(2 * M + 1);c1++)
    resul = resul + A[1][c1] * A[1][c1];

  resul = (double)sqrt((double)resul);
  return(resul);
 }

void Perfaz_Soma1(double ***A, double ***B, long N, long M)
 {
  long ca,cb;

  for (ca=1;ca<(N + (long)1);ca++)
    for (cb=1;cb<(M + (long)1);cb++)
      A[1][ca][cb] = A[1][ca][cb] + B[1][ca][cb];
 }

void Perfaz_Soma2(double **a, double **b, long M)
 {
  long ca;
  for (ca=1;ca<(2 * M + (long)1);ca++)
    a[1][ca] = a[1][ca] + b[1][ca];
 }

void Gerokz(double ***So,double **so,long N,long M,double Dx,double Dy)
 {
  extern double interv_x,interv_y;
  extern double rho,zzero,delta;
  extern long td_x,td_y;

  double u,v;
  double Dv,Du,fator;
  long i,j;

  Du = ((double)2. * PI / (((double)N)*Dx));
  Dv = ((double)2. * PI / (((double)M)*Dy));

  for(j=1;j<((M/2) +1);j++)
   {
    v = (double)(j - 1) * Dv;
    for(i=1;i<((N/2) +1);i++)
     {
      u = (double)(i - 1) * Du;
      fator = Fator(u,v,N,M);
      So[1][i][2*j -1] = fator*So[1][i][2*j -1];
      So[1][i][2*j]    = fator*So[1][i][2*j];
     }
    for(i=((N/2) +1);i<(N +1);i++)
     {
      u = (double)(i - N - 1) * Du;
      fator = Fator(u,v,N,M);
      So[1][i][2*j -1] = fator*So[1][i][2*j -1];
      So[1][i][2*j]    = fator*So[1][i][2*j];
     }
   }

  v = (double)(M/2) * Dv;
    for(i=1;i<((N/2) +1);i++)
     {
      u = (double)(i - 1) * Du;
      fator = Fator(u,v,N,M);
      so[1][2 * i - 1] = fator*so[1][2 * i - 1];
      so[1][2 * i] = fator*so[1][2 * i];
     }
    for(i=((N/2) +1);i<(N +1);i++)
     {
      u = (double)(i - N - 1) * Du;
      fator = Fator(u,v,N,M);
      so[1][2 * i - 1] = fator*so[1][2 * i - 1];
      so[1][2 * i] = fator*so[1][2 * i];
     }
 }

double Fator(double U, double V, long N, long M)
 {
  extern double interv_x,interv_y;
  extern double rho,zzero,delta;
  extern long td_x,td_y;

  double k,resul,esc;

  k = (double)sqrt((double)(U*U + V*V));
  esc = (double)2. / ((double)N * (double)M);

  resul = esc * (double)(-2.) * PI * G * rho * (double)exp((double)((double)(-1.) * k * zzero));

  return(resul);
 }

double maxover(double ***So,double **so,long N,long M,double Dx,double Dy)
 {
  extern double interv_x,interv_y;
  extern double rho,zzero,delta;
  extern long td_x,td_y;

  double u,v;
  double Dv,Du,fator,aux;
  long i,j;
  double cand_max,Max=(double)0.;

  Du = ((double)2. * PI / (((double)N)*Dx));
  Dv = ((double)2. * PI / (((double)M)*Dy));

  for(j=1;j<((M/2) +1);j++)
   {
    v = (double)(j - 1) * Dv;
    for(i=1;i<((N/2) +1);i++)
     {
      u = (double)(i - 1) * Du;
      fator = FatorC(u,v);
      aux = fator * fator;
      cand_max = aux * (So[1][i][2*j -1] * So[1][i][2*j -1] + So[1][i][2*j] * So[1][i][2*j]);
      if (cand_max > Max) Max = cand_max;
     }
    for(i=((N/2) +1);i<(N +1);i++)
     {
      u = (double)(i - N - 1) * Du;
      fator = FatorC(u,v);
      aux = fator * fator;
      cand_max = aux * (So[1][i][2*j -1] * So[1][i][2*j -1] + So[1][i][2*j] * So[1][i][2*j]);
      if (cand_max > Max) Max = cand_max;
     }
   }

  v = (double)(M/2) * Dv;
    for(i=1;i<((N/2) +1);i++)
     {
      u = (double)(i - 1) * Du;
      fator = FatorC(u,v);
      aux = fator * fator;
      cand_max = aux * (so[1][2*i -1] * so[1][2*i -1] + so[1][2*i] * so[1][2*i]);
      if (cand_max > Max) Max = cand_max;
     }
    for(i=((N/2) +1);i<(N +1);i++)
     {
      u = (double)(i - N - 1) * Du;
      fator = FatorC(u,v);
      aux = fator * fator;
      cand_max = aux * (so[1][2*i -1] * so[1][2*i -1] + so[1][2*i] * so[1][2*i]);
      if (cand_max > Max) Max = cand_max;
     }

   Max = (double)sqrt((double)Max);
   return(Max);
 }


double FatorC(double U, double V)
 {
  extern double interv_x,interv_y;
  extern double rho,zzero,delta;
  extern long td_x,td_y;

  double k,resul,esc;

  k = (double)sqrt((double)(U*U + V*V));

  resul = (double)exp((double)((double)(-1.) * k * zzero));

  return(resul);
 }

