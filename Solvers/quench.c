/*****************************************************************************
  Quench is a library for calculation of quenching data.
  See quench.h for details.
  Last updated: 4 february 1995 , by KN
*****************************************************************************/
#include "quench.h"
#include <math.h>
#include <stdio.h>

double arctan_local(double complex z) {
  #define PI 3.14159265359
  double angle;

  if (creal(z) != 0.0) {
    angle= carg(z); 
    /*
    if (creal(z)<0.0) 
      angle += PI;
      */
    if (angle >PI) 
      angle -= 2.0*PI;
   } else {
    if (cimag(z)>0.0) {
      angle= PI/2.0;
    } else {
      angle= PI/(-2.0);
    } 
   }
  return (360.0/2.0/PI*angle);
} /* arctan_local */

void compamppha(int n, double **P, double *amp, double *phase) {

   int   i;
   double complex z;

   for(i=0;i<n;i++) {
     z= P[i][0] + P[i][1]*I;
     amp[i]= cabs(z);

/* return (360.0/2.0/3.14159265359*atan(h,g)) */
/*     phase[i]= (57.2957795131*carg(z)); */
   phase[i]= arctan_local(z);
   }
} /* compamppha */

void stopdata(int n, int m, double **IP, double *x,
	      double *q, double *fi, double *qd, double *fid) {

   int     i;
   double  gd, hd;
   double complex z;

   for(i=0;i<n;i++) {
     z= -IP[0][i] - IP[1][i]*I;
     q[i]=  1.0/cabs(z);
     fi[i]= (57.2957795131*carg(z));
   }

   gd= 0.0; hd= 0.0;
   for(i=0;i<n;i++) {
     gd += IP[0][i]*x[i];
     hd += IP[1][i]*x[i];
   }

   z= gd + hd*I;
   *qd=  -x[m-1]/cabs(z);
   *fid= (57.2957795131*carg(z));

} /* stopdata */

