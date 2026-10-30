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
  case 0: return ((-(y)));
break;
  case 1: return (x);
break;
  } /* switch */
} /* reac */
void eval(double S[equa], double F[equa]) {
 double x;
 double y;
  x=S[0];
  y=S[1];
  F[0] = (-(y));
  F[1] = x;
} /* eval */
void calc_jac(double S[equa], double tk, double R[equa], double Jac[equa+1][equa+1]) 
{
 double x;
 double y;
  x=S[0];
  y=S[1];
Jac[1][1]=tk*(-0.000000e+00)-(1+4*R[0]);
Jac[1][2]=tk*(-1.000000e+00);
Jac[2][1]=tk*(1.000000e+00);
Jac[2][2]=tk*(0.000000e+00)-(1+4*R[1]);
}
double reac_F(double S[equa], double K[equa], double R[equa], double tk, int l) {
 double x;
 double y;
  x=S[0];
  y=S[1];
switch (l) {
case 1:
return (tk*((-(y)))-S[0]*(1+4*R[0])+K[0]);
break;
case 2:
return (tk*(x)-S[1]*(1+4*R[1])+K[1]);
break;
}
}
void calc_jac2(double S[equa], double Jac[equa][equa]) 
{
double x;
double y;
x=S[0];
y=S[1];
Jac[0][0]=-0.000000e+00;
Jac[0][1]=-1.000000e+00;
Jac[1][0]=1.000000e+00;
Jac[1][1]=0.000000e+00;
}
