/* ########################
*/

#include <math.h>
#include "model.h"
void init_diff_const(void) {
} /* init_diff_const */

double reac(double* S, int l) {
 double A;
 double P;
  A=S[0];
  P=S[1];
  switch (l) {
  case 0: return ((-1.000000e+00*A));
break;
  case 1: return (A);
break;
  } /* switch */
} /* reac */
void eval(double S[equa], double F[equa]) {
 double A;
 double P;
  if (S[0]<0.0) S[0]=0.0;
  A=S[0];
  if (S[1]<0.0) S[1]=0.0;
  P=S[1];
  F[0] = (-1.000000e+00*A);
  F[1] = A;
} /* eval */
void calc_jac(double S[equa], double tk, double R[equa], double Jac[equa+1][equa+1]) 
{
 double A;
 double P;
  if (S[0]<0.0) S[0]=0.0;
  A=S[0];
  if (S[1]<0.0) S[1]=0.0;
  P=S[1];
Jac[1][1]=tk*(-1.000000e+00)-(1+4*R[0]);
Jac[1][2]=tk*(0.000000e+00);
Jac[2][1]=tk*(1.000000e+00);
Jac[2][2]=tk*(0.000000e+00)-(1+4*R[1]);
}
double reac_F(double S[equa], double K[equa], double R[equa], double tk, int l) {
 double A;
 double P;
  A=S[0];
  P=S[1];
switch (l) {
case 1:
return (tk*((-1.000000e+00*A))-S[0]*(1+4*R[0])+K[0]);
break;
case 2:
return (tk*(A)-S[1]*(1+4*R[1])+K[1]);
break;
}
}
void calc_jac2(double S[equa], double Jac[equa][equa]) 
{
double A;
double P;
A=S[0];
P=S[1];
Jac[0][0]=-1.000000e+00;
Jac[0][1]=0.000000e+00;
Jac[1][0]=1.000000e+00;
Jac[1][1]=0.000000e+00;
}
