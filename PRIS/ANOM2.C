/* anom2 : Rotina para o calculo de anomalia gravimetrica
	   devida a um prisma vertical apresentando contraste de
	   densidade variando quadraticamente com a
	   profundidade. (EQUACAO APROXIMADA)

*/

#include<stdio.h>
#include<math.h>

double anom2(double x0,double y0,double xp,double yp,double z1,double z2)
{
 extern double interv_x,interv_y,a0,a1,a2;

 double a;
 double gamma=(double)6.670e-3;
 double gpris;

 double x1,xx1,y1,yy1,r1,r2,p1,p2,p3;

 a = interv_x * interv_y;
 xx1 = fabs(x0 - xp);
 yy1 = fabs(y0 - yp);

 x1 = xx1;
 y1 = yy1;

 if(x1 == (double)0.) x1 = (double)0.0000001;
 if(y1 == (double)0.) y1 = (double)0.0000001;

 r1 = sqrt(x1*x1 + y1*y1 + z1*z1);
 r2 = sqrt(x1*x1 + y1*y1 + z2*z2);

 if(r1 == (double)0.) r1 = (double)0.0000001;
 if(r2 == (double)0.) r2 = (double)0.0000001;

 p1 = a0 * ( ( (double)1.0 / r1 ) - ( (double)1.0 / r2) );
 p2 = a1 * ( z1/r1 - z2/r2 + log((r2+z2) / (r1+z1)) );
 p3 = a2 * ( ((double)2.0 * (r2 - r1)) + (z1 * z1 / r1) - (z2 * z2 / r2) );

 gpris = gamma * a * ( p1 + p2 + p3 );

 return(gpris);

}
