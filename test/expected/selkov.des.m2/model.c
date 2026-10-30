/* ########################
*/

#include <math.h>
#include "model.h"
void init_diff_const(void) {
} /* init_diff_const */

double reac(double* S, int l) {
 double x;
 double y;
  x=S[1];
  y=S[2];
  switch (l) {
  case 0: return ((1.000000e+00-(x*(y*y*y))));
break;
  case 1: return ((3.000000e-01*((x*(y*y*y))-y)));
break;
  } /* switch */
} /* reac */
void eval(double S[equa], double F[equa]) {
 double x;
 double y;
  x=S[0];
  y=S[1];
  F[0] = (1.000000e+00-(x*(y*y*y)));
  F[1] = (3.000000e-01*((x*(y*y*y))-y));
} /* eval */
void calc_jac(double S[equa], double tk, double R[equa], double Jac[equa+1][equa+1]) 
{
 double x;
 double y;
  x=S[0];
  y=S[1];
Jac[1][1]=tk*((-((y*y*y))))-(1+4*R[0]);
Jac[1][2]=tk*((-((x*((y*y)*3.000000e+00)))));
Jac[2][1]=tk*((3.000000e-01*(y*y*y)));
Jac[2][2]=tk*((3.000000e-01*((x*((y*y)*3.000000e+00))-1.000000e+00
)))-(1+4*R[1]);
}
double reac_F(double S[equa], double K[equa], double R[equa], double tk, int l) {
 double x;
 double y;
  x=S[0];
  y=S[1];
switch (l) {
case 1:
return (tk*((1.000000e+00-(x*(y*y*y))))-S[0]*(1+4*R[0])+K[0]);
break;
case 2:
return (tk*((3.000000e-01*((x*(y*y*y))-y)))-S[1]*(1+4*R[1])+K[1]);
break;
}
}
void calc_jac2(double S[equa], double Jac[equa][equa]) 
{
double x;
double y;
x=S[0];
y=S[1];
Jac[0][0]=(-((y*y*y)));
Jac[0][1]=(-((x*((y*y)*3.000000e+00))));
Jac[1][0]=(3.000000e-01*(y*y*y));
Jac[1][1]=(3.000000e-01*((x*((y*y)*3.000000e+00))-1.000000e+00
));
}
