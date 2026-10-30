/* ########################
*/

#include <math.h>
#include "model.h"
void init_diff_const(void) {
} /* init_diff_const */

double reac(double* S, int l) {
 double x;
 double t;
  x=S[1];
  t=S[2];
  switch (l) {
  case 0: return (cos(t));
break;
  case 1: return (1.000000e+00);
break;
  } /* switch */
} /* reac */
void eval(double S[equa], double F[equa]) {
 double x;
 double t;
  x=S[0];
  t=S[1];
  F[0] = cos(t);
  F[1] = 1.000000e+00;
} /* eval */
void calc_jac(double S[equa], double tk, double R[equa], double Jac[equa+1][equa+1]) 
{
 double x;
 double t;
  x=S[0];
  t=S[1];
Jac[1][1]=tk*(0.000000e+00)-(1+4*R[0]);
Jac[1][2]=tk*((-(sin(t))));
Jac[2][1]=tk*(0.000000e+00);
Jac[2][2]=tk*(0.000000e+00)-(1+4*R[1]);
}
double reac_F(double S[equa], double K[equa], double R[equa], double tk, int l) {
 double x;
 double t;
  x=S[0];
  t=S[1];
switch (l) {
case 1:
return (tk*(cos(t))-S[0]*(1+4*R[0])+K[0]);
break;
case 2:
return (tk*(1.000000e+00)-S[1]*(1+4*R[1])+K[1]);
break;
}
}
void calc_jac2(double S[equa], double Jac[equa][equa]) 
{
double x;
double t;
x=S[0];
t=S[1];
Jac[0][0]=0.000000e+00;
Jac[0][1]=(-(sin(t)));
Jac[1][0]=0.000000e+00;
Jac[1][1]=0.000000e+00;
}
