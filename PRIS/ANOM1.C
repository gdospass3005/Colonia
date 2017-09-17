/* anom1 : Rotina para o calculo de anomalia gravimetrica
	   devida a um prisma vertical apresentando contraste de
	   densidade variando quadraticamente com a
	   profundidade. (EQUACAO EXATA)    

*/

#include<stdio.h>
#include<math.h>

double anom1(double x0,double y0,double xp,double yp,double z1,double z2)
{
 extern double interv_x,interv_y,a0,a1,a2;

 double xx1,yy1;
 double t,w;
 double x1,x2,y1,y2;
 double r1,r2,r3,r4,r5,r6,r7,r8;
 double f11,f12,f13,f14,f21,f22,f23,f24,f31,f32,f33,f34;
 double f41,f42,f43,f44,f51,f52,f53,f54,f61,f62,f63,f64,f71,f72,f73,f74;
 double f81,f82,f83,f84,f91,f92,f93,f94,f101,f102,f103,f104;
 double h71,h72,h81,h82,h91,h92,h101,h102;
 double b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14;
 double gamma=(double)6.670e-3;
 double g1=(double)0.,g2=(double)0.,g3=(double)0.;
 double gpris;

 xx1 = fabs(x0 - xp);
 yy1 = fabs(y0 - yp);

 t = ( interv_x / (double)2. );
 w = ( interv_y / (double)2. );

 x1 = xx1 + t;
 x2 = xx1 - t;
 y1 = yy1 + w;
 y2 = yy1 - w;

 if(x1 == (double)0.) x1 = (double)0.0000001;
 if(x2 == (double)0.) x2 = (double)0.0000001;
 if(y1 == (double)0.) y1 = (double)0.0000001;
 if(y2 == (double)0.) y2 = (double)0.0000001;
    
 r1 = sqrt( x2 * x2 + y2 * y2 + z1 * z1 ) ;
 r2 = sqrt( x2 * x2 + y2 * y2 + z2 * z2 ) ;
 r3 = sqrt( x2 * x2 + y1 * y1 + z1 * z1 ) ;
 r4 = sqrt( x2 * x2 + y1 * y1 + z2 * z2 ) ;
 r5 = sqrt( x1 * x1 + y2 * y2 + z1 * z1 ) ;
 r6 = sqrt( x1 * x1 + y2 * y2 + z2 * z2 ) ;
 r7 = sqrt( x1 * x1 + y1 * y1 + z1 * z1 ) ;
 r8 = sqrt( x1 * x1 + y1 * y1 + z2 * z2 ) ;

 f11 = x2 * y2 / (z2 * r2) ;
 f12 = x2 * y1 / (z2 * r4) ;    
 f13 = x1 * y2 / (z2 * r6) ;    
 f14 = x1 * y1 / (z2 * r8) ;    
  b1 = atan(f11) - atan(f12) - atan(f13) + atan(f14) ; 
    
 f21 = x2 * y2 / (z1 * r1) ;
 f22 = x2 * y1 / (z1 * r3) ;    
 f23 = x1 * y2 / (z1 * r5) ;    
 f24 = x1 * y1 / (z1 * r7) ;    
  b2 = atan(f21) - atan(f22) - atan(f23) + atan(f24) ;
    
 f31 = (r2 - y2) / (r2 + y2) ;
 f32 = (r1 - y2) / (r1 + y2) ; 
 f33 = (r4 - y1) / (r4 + y1) ;
 f34 = (r3 - y1) / (r3 + y1) ;
  b3 = log(f31 * f34 / (f32 * f33)) ;

 f41 = (r5 - y2) / (r5 + y2) ;
 f42 = (r6 - y2) / (r6 + y2) ; 
 f43 = (r7 - y1) / (r7 + y1) ;
 f44 = (r8 - y1) / (r8 + y1) ;
  b4 = log(f41 * f44 / (f42 * f43)) ;

 f51 = (r2 - x2) / (r2 + x2) ;
 f52 = (r1 - x2) / (r1 + x2) ; 
 f53 = (r6 - x1) / (r6 + x1) ;
 f54 = (r5 - x1) / (r5 + x1) ;
  b5 = log(f51 * f54 / (f52 * f53)) ;

 f61 = (r3 - x2) / (r3 + x2) ;
 f62 = (r4 - x2) / (r4 + x2) ; 
 f63 = (r7 - x1) / (r7 + x1) ;
 f64 = (r8 - x1) / (r8 + x1) ;
  b6 = log(f61 * f64 / (f62 * f63)) ;

 if ( a0 != (double)0.)
  g1 = gamma * a0 * ( (z2 * b1) - (z1 * b2) + ((x2 / (double)2.) * b3) +
   ((x1 / (double)2.) * b4) +
   ((y2 / (double)2.) * b5) + ((y1 / (double)2.) * b6) ) ;

 if ( a1 != (double)0. || a2 != (double)0. )
  {

   f71 = y2 * z2 / (r2 * x2) ;
   f72 = y2 * z1 / (r1 * x2) ;
   f73 = y1 * z2 / (r4 * x2) ;
   f74 = y1 * z1 / (r3 * x2) ;
   h71 = (f71 - f72) / ( (double)1. + (f71 * f72) ) ;
   h72 = (f73 - f74) / ( (double)1. + (f73 * f74) ) ;
    b7 = atan( (h71 - h72) / ( (double)1. + (h71 * h72) ) ) ;

   f81 = y2 * z2 / (r6 * x1) ;
   f82 = y2 * z1 / (r5 * x1) ;
   f83 = y1 * z2 / (r8 * x1) ;
   f84 = y1 * z1 / (r7 * x1) ;
   h81 = (f81 - f82) / ( (double)1. + (f81 * f82) ) ;
   h82 = (f83 - f84) / ( (double)1. + (f83 * f84) ) ;
    b8 = atan( (h81 - h82) / ( (double)1. + (h81 * h82) ) ) ;

   f91 = x2 * z2 / (r2 * y2) ;
   f92 = x2 * z1 / (r1 * y2) ;
   f93 = x1 * z2 / (r6 * y2) ;
   f94 = x1 * z1 / (r5 * y2) ;
   h91 = (f91 - f92) / ( (double)1. + (f91 * f92) ) ;
   h92 = (f93 - f94) / ( (double)1. + (f93 * f94) ) ;
    b9 = atan( (h91 - h92) / ( (double)1. + (h91 * h92) ) ) ;

   f101 = x2 * z2 / (r4 * y1) ;
   f102 = x2 * z1 / (r3 * y1) ;
   f103 = x1 * z2 / (r8 * y1) ;
   f104 = x1 * z1 / (r7 * y1) ;
   h101 = (f101 - f102) / ( (double)1. + (f101 * f102) ) ;
   h102 = (f103 - f104) / ( (double)1. + (f103 * f104) ) ;
    b10 = atan( (h101 - h102) / ( (double)1. + (h101 * h102) ) ) ;

   b11 = log( (r2 + z2) / (r1 + z1) ) ;
   b12 = log( (r3 + z1) / (r4 + z2) ) ;
   b13 = log( (r5 + z1) / (r6 + z2) ) ;
   b14 = log( (r8 + z2) / (r7 + z1) ) ;

   if (a1 != (double)0.)
    g2 = (gamma /(double)2.) * a1 * ( (z2 * z2 * b1) - (z1 * z1 * b2) - 
     (x2 * x2 *b7) + (x1 * x1 * b8) - (y2 * y2 * b9) + (y1 * y1 * b10) +
     ( (double)2. * ( (x2 * y2 *b11) + (x2 * y1 * b12) + (x1 * y2 * b13) +
     (x1 * y1 * b14) ) ) ) ;
 
   if (a2 != (double)0.)
    g3 = (gamma /(double)3.) * a2 * ( (z2 * z2 * z2 * b1) - 
     (z1 * z1 * z1 * b2) - ( ((x2 * x2 * x2) / (double)2.) * b3 ) -
     ( ((x1 * x1 * x1) / (double)2.) * b4 ) -
     ( ((y2 * y2 * y2) / (double)2.) * b5 ) -
     ( ((y1 * y1 * y1) / (double)2.) * b6 ) +
     (double)2. * ( (x2 * y2 * (r2 - r1)) + (x2 * y1 * (r3 - r4)) +
     (x1 * y2 * (r5 - r6)) + (x1 * y1 * (r8 - r7)) ) ) ;

  }

 gpris = g1 + g2 + g3;

 return(gpris);

}
