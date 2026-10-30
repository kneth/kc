/* ########################
*/

#include <math.h>
#include "model.h"
void init_diff_const(void) {
} /* init_diff_const */

double reac(double* S, int l) {
 double A;
 double X;
 double Y;
  A=S[0];
  X=S[1];
  Y=S[2];
  switch (l) {
  case 0: return (((-1.000000e+00*((1.000000e-01*A)*X))+((1.000000e-01
*Y)*A)));
break;
  case 1: return ((((1.000000e-01*A)*X)+(-1.000000e+00*((1.010000e-01
*X)*Y))));
break;
  case 2: return ((((1.010000e-01*X)*Y)+(-1.000000e+00*((1.000000e-01
*Y)*A))));
break;
  } /* switch */
} /* reac */
void eval(double S[equa], double F[equa]) {
 double A;
 double X;
 double Y;
  if (S[0]<0.0) S[0]=0.0;
  A=S[0];
  if (S[1]<0.0) S[1]=0.0;
  X=S[1];
  if (S[2]<0.0) S[2]=0.0;
  Y=S[2];
  F[0] = ((-1.000000e+00*((1.000000e-01*A)*X))+((1.000000e-01
*Y)*A));
  F[1] = (((1.000000e-01*A)*X)+(-1.000000e+00*((1.010000e-01
*X)*Y)));
  F[2] = (((1.010000e-01*X)*Y)+(-1.000000e+00*((1.000000e-01
*Y)*A)));
} /* eval */
void calc_jac(double S[equa], double tk, double R[equa], double Jac[equa+1][equa+1]) 
{
 double A;
 double X;
 double Y;
  if (S[0]<0.0) S[0]=0.0;
  A=S[0];
  if (S[1]<0.0) S[1]=0.0;
  X=S[1];
  if (S[2]<0.0) S[2]=0.0;
  Y=S[2];
Jac[1][1]=tk*(((-1.000000e+00*(1.000000e-01*X))+(1.000000e-01
*Y)))-(1+4*R[0]);
Jac[1][2]=tk*((-1.000000e+00*(1.000000e-01*A)));
Jac[1][3]=tk*((1.000000e-01*A));
Jac[2][1]=tk*((1.000000e-01*X));
Jac[2][2]=tk*(((1.000000e-01*A)+(-1.000000e+00*(1.010000e-01
*Y))))-(1+4*R[1]);
Jac[2][3]=tk*((-1.000000e+00*(1.010000e-01*X)));
Jac[3][1]=tk*((-1.000000e+00*(1.000000e-01*Y)));
Jac[3][2]=tk*((1.010000e-01*Y));
Jac[3][3]=tk*(((1.010000e-01*X)+(-1.000000e+00*(1.000000e-01
*A))))-(1+4*R[2]);
}
reac_F(double S[equa], double K[equa], double R[equa], double tk, int l) {
 double A;
 double X;
 double Y;
  A=S[0];
  X=S[1];
  Y=S[2];
switch (l) {
case 1:
return (tk*(((-1.000000e+00*((1.000000e-01*A)*X))+((1.000000e-01
*Y)*A)))-S[0]*(1+4*R[0])+K[0]);
break;
case 2:
return (tk*((((1.000000e-01*A)*X)+(-1.000000e+00*((1.010000e-01
*X)*Y))))-S[1]*(1+4*R[1])+K[1]);
break;
case 3:
return (tk*((((1.010000e-01*X)*Y)+(-1.000000e+00*((1.000000e-01
*Y)*A))))-S[2]*(1+4*R[2])+K[2]);
break;
}
}
void calc_jac2(double S[equa], double Jac[equa][equa]) 
{
double A;
double X;
double Y;
A=S[0];
X=S[1];
Y=S[2];
Jac[0][0]=((-1.000000e+00*(1.000000e-01*X))+(1.000000e-01
*Y));
Jac[0][1]=(-1.000000e+00*(1.000000e-01*A));
Jac[0][2]=(1.000000e-01*A);
Jac[1][0]=(1.000000e-01*X);
Jac[1][1]=((1.000000e-01*A)+(-1.000000e+00*(1.010000e-01
*Y)));
Jac[1][2]=(-1.000000e+00*(1.010000e-01*X));
Jac[2][0]=(-1.000000e+00*(1.000000e-01*Y));
Jac[2][1]=(1.010000e-01*Y);
Jac[2][2]=((1.010000e-01*X)+(-1.000000e+00*(1.000000e-01
*A)));
}
