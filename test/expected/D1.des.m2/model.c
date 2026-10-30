/* ########################
*/

#include <math.h>
#include "model.h"
void init_diff_const(void) {
} /* init_diff_const */

double reac(double* S, int l) {
 double y1;
 double y2;
 double y3;
  y1=S[1];
  y2=S[2];
  y3=S[3];
  switch (l) {
  case 0: return ((2.000000e-01*(y2-y1)));
break;
  case 1: return ((((1.000000e+01*y1)-((6.000000e+01-(1.250000e-01
*y3))*y2))+(1.250000e-01*y3)));
break;
  case 2: return (1.000000e+00);
break;
  } /* switch */
} /* reac */
void eval(double S[equa], double F[equa]) {
 double y1;
 double y2;
 double y3;
  y1=S[0];
  y2=S[1];
  y3=S[2];
  F[0] = (2.000000e-01*(y2-y1));
  F[1] = (((1.000000e+01*y1)-((6.000000e+01-(1.250000e-01
*y3))*y2))+(1.250000e-01*y3));
  F[2] = 1.000000e+00;
} /* eval */
void calc_jac(double S[equa], double tk, double R[equa], double Jac[equa+1][equa+1]) 
{
 double y1;
 double y2;
 double y3;
  y1=S[0];
  y2=S[1];
  y3=S[2];
Jac[1][1]=tk*(-2.000000e-01)-(1+4*R[0]);
Jac[1][2]=tk*(2.000000e-01);
Jac[1][3]=tk*(0.000000e+00);
Jac[2][1]=tk*(1.000000e+01);
Jac[2][2]=tk*((-((6.000000e+01-(1.250000e-01*y3)))))-(1+4*R[1]);
Jac[2][3]=tk*(((-((-1.250000e-01*y2)))+1.250000e-01));
Jac[3][1]=tk*(0.000000e+00);
Jac[3][2]=tk*(0.000000e+00);
Jac[3][3]=tk*(0.000000e+00)-(1+4*R[2]);
}
reac_F(double S[equa], double K[equa], double R[equa], double tk, int l) {
 double y1;
 double y2;
 double y3;
  y1=S[0];
  y2=S[1];
  y3=S[2];
switch (l) {
case 1:
return (tk*((2.000000e-01*(y2-y1)))-S[0]*(1+4*R[0])+K[0]);
break;
case 2:
return (tk*((((1.000000e+01*y1)-((6.000000e+01-(1.250000e-01
*y3))*y2))+(1.250000e-01*y3)))-S[1]*(1+4*R[1])+K[1]);
break;
case 3:
return (tk*(1.000000e+00)-S[2]*(1+4*R[2])+K[2]);
break;
}
}
void calc_jac2(double S[equa], double Jac[equa][equa]) 
{
double y1;
double y2;
double y3;
y1=S[0];
y2=S[1];
y3=S[2];
Jac[0][0]=-2.000000e-01;
Jac[0][1]=2.000000e-01;
Jac[0][2]=0.000000e+00;
Jac[1][0]=1.000000e+01;
Jac[1][1]=(-((6.000000e+01-(1.250000e-01*y3))));
Jac[1][2]=((-((-1.250000e-01*y2)))+1.250000e-01);
Jac[2][0]=0.000000e+00;
Jac[2][1]=0.000000e+00;
Jac[2][2]=0.000000e+00;
}
