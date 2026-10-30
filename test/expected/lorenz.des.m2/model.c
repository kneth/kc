/* ########################
*/

#include <math.h>
#include "model.h"
void init_diff_const(void) {
} /* init_diff_const */

double reac(double* S, int l) {
 double X;
 double Y;
 double Z;
  X=S[1];
  Y=S[2];
  Z=S[3];
  switch (l) {
  case 0: return (((-1.600000e+01*X)+(1.600000e+01*Y)));
break;
  case 1: return ((((4.592000e+01*X)-(X*Z))-Y));
break;
  case 2: return (((X*Y)-(4.000000e+00*Z)));
break;
  } /* switch */
} /* reac */
void eval(double S[equa], double F[equa]) {
 double X;
 double Y;
 double Z;
  X=S[0];
  Y=S[1];
  Z=S[2];
  F[0] = ((-1.600000e+01*X)+(1.600000e+01*Y));
  F[1] = (((4.592000e+01*X)-(X*Z))-Y);
  F[2] = ((X*Y)-(4.000000e+00*Z));
} /* eval */
void calc_jac(double S[equa], double tk, double R[equa], double Jac[equa+1][equa+1]) 
{
 double X;
 double Y;
 double Z;
  X=S[0];
  Y=S[1];
  Z=S[2];
Jac[1][1]=tk*(-1.600000e+01)-(1+4*R[0]);
Jac[1][2]=tk*(1.600000e+01);
Jac[1][3]=tk*(0.000000e+00);
Jac[2][1]=tk*((4.592000e+01-Z));
Jac[2][2]=tk*(-1.000000e+00)-(1+4*R[1]);
Jac[2][3]=tk*((-(X)));
Jac[3][1]=tk*(Y);
Jac[3][2]=tk*(X);
Jac[3][3]=tk*(-4.000000e+00)-(1+4*R[2]);
}
reac_F(double S[equa], double K[equa], double R[equa], double tk, int l) {
 double X;
 double Y;
 double Z;
  X=S[0];
  Y=S[1];
  Z=S[2];
switch (l) {
case 1:
return (tk*(((-1.600000e+01*X)+(1.600000e+01*Y)))-S[0]*(1+4*R[0])+K[0]);
break;
case 2:
return (tk*((((4.592000e+01*X)-(X*Z))-Y))-S[1]*(1+4*R[1])+K[1]);
break;
case 3:
return (tk*(((X*Y)-(4.000000e+00*Z)))-S[2]*(1+4*R[2])+K[2]);
break;
}
}
void calc_jac2(double S[equa], double Jac[equa][equa]) 
{
double X;
double Y;
double Z;
X=S[0];
Y=S[1];
Z=S[2];
Jac[0][0]=-1.600000e+01;
Jac[0][1]=1.600000e+01;
Jac[0][2]=0.000000e+00;
Jac[1][0]=(4.592000e+01-Z);
Jac[1][1]=-1.000000e+00;
Jac[1][2]=(-(X));
Jac[2][0]=Y;
Jac[2][1]=X;
Jac[2][2]=-4.000000e+00;
}
