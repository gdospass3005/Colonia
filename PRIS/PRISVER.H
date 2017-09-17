double gbox(double x0, double y0, double z0, double x1, double y1, double z1,
double x2, double y2, double z2, double rho);
double anom1(double x0,double y0,double x1,double y1,double z1,double z2);
double anom2(double x0,double y0,double x1,double y1,double z1,double z2);

void load_h(char *arq);
void load_g(char *arq);
void load_coef(char *arq);
void save(char *arq2, double **G);
void saverms(char *arq2, double *RMSH, double *RMSG, double *MAXERRG, double *MAXERRH);
void avalia(double **H, double **G, long iter);
void flat(double **G, double **H);
void loop(double **G, double **H);
void nrerror(char error_text[]);
double *dvector(long nl, long nh);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
void free_dvector(double *v, long nl, long nh);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
