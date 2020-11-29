/*****************************************************************************
  Quench is a library for calculation of quenching data.
  See quench.h for details.
  Last updated: 4 february 1995 , by KN
*****************************************************************************/
#include "complex.h"
#include "matrix.h"
#include <math.h>
#include <stdio.h>

double arctan_local(Complex z) {
  #define PI 3.14159265359
  double angle;

  if (z.re != 0.0) {
    angle= ComplexArg(z); 
    /*
    if (z.re<0.0) 
      angle += PI;
      */
    if (angle >PI) 
      angle -= (2.0*PI);
   } else {
    if (z.im>0.0) {
      angle= PI/2.0;
    } else {
      angle= PI/(-2.0);
    } 
   }
  return (360.0/2.0/PI*angle);
} /* arctan_local */

void compamppha(int n, double **P, double *amp, double *phase) {

   int   i;
   Complex z;
   double hlp;

   for(i=0;i<n;i++) {
     ComplexAssign(P[i][0],P[i][1],&z);
     amp[i]= ComplexNorm(z);

/* return (360.0/2.0/3.14159265359*atan(h,g)) */
/*     phase[i]= (57.2957795131*ComplexArg(z)); */
   phase[i]= arctan_local(z);
   }
} /* compamppha */

void stopdata(int n, int m, double **IP, double *x,
	      double *q, double *fi, double *qd, double *fid) {

   int     i;
   double  gd, hd;
   Complex z;

   for(i=0;i<n;i++) {
     ComplexAssign(-IP[0][i],-IP[1][i],&z);
     q[i]=  1.0/ComplexNorm(z);
     fi[i]= (57.2957795131*ComplexArg(z));
   }

   gd= 0.0; hd= 0.0;
   for(i=0;i<n;i++) {
     gd += (IP[0][i]*x[i]);
     hd += (IP[1][i]*x[i]);
   }

   ComplexAssign(gd,hd,&z);
   *qd=  -x[m-1]/ComplexNorm(z);
   *fid= (57.2957795131*ComplexArg(z));

} /* stopdata */

