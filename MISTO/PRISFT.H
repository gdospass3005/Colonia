#ifndef _NR_H_
#define _NR_H_
#ifndef _NR_UTILS_H_
#define _NR_UTILS_H_

#include <stdio.h>

double gbox(double x0, double y0, double z0, double x1, double y1, double z1,
double x2, double y2, double z2, double rho);
double anom1(double x0,double y0,double x1,double y1,double z1,double z2);
double anom2(double x0,double y0,double x1,double y1,double z1,double z2);

double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void free_d3tensor(double ***t, long nrl, long nrh, long ncl, long nch,
        long ndl, long ndh);

void load_h(char *arq);
void load_g(char *arq);
void load_coef_ift(char *arq);
void save(char *arq2, double **G);
void saverms(char *arq2, double *RMSH, double *RMSG, double *MAXERRG, double *MAXERRH);
void avalia(double **H, double **G, long iter);
void avaliaft(double **H, double **G);
void flat(double **G, double **H);
void loop(double **G, double **H);
void loopft(double **G, double **H);
void mudaoriz(double **,double);

void nrerror(char error_text[]);
double *dvector(long nl, long nh);
double **dmatrix(long nrl, long nrh, long ncl, long nch);


void Perfaz_Subt1(double ***,double ***,long,long);
void Perfaz_Subt2(double **, double **, long);
void parcel(double ***,double **,long,long,double,double,long,long);
double ffator(double, double, long, long);
void Perfaz_Soma1(double ***,double ***,long,long);
void Perfaz_Soma2(double **, double **, long);
void Gerokz(double ***,double **,long,long,double,double);
void Gerokzi(double ***,double **,long,long,double,double);
double Fator2(double, double, long, long);
double Fator(double, double, long, long);
double FatorC(double, double);
double fator(double, double);
double maxover(double ***,double **,long,long,double,double);
void Filtro(double ***,double **,long,long,double,double,double,double);
double filtro(double, double, double, double);

double pythag(double a, double b);
void svbksb(double **u, double w[], double **v, int m, int n, double b[],
        double x[]);
void svdcmp(double **a, int m, int n, double w[], double **v);
void svdfit(double x[], double y[], double sig[], int ndata, double a[],
        int ma, double **u, double **v, double w[], double *chisq,
        void (*funcs)(double, double [], int));
void svdvar(double **v, int ma, double w[], double **cvm);
void fourn(double data[], unsigned long nn[], int ndim, int isign);
void rlft3(double ***data, double **speq, unsigned long nn1,
        unsigned long nn2, unsigned long nn3, int isign);



#endif /* _NR_UTILS_H_ */
#endif /* _NR_H_ */


