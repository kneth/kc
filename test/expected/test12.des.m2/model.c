/* ########################
*/

#include <math.h>
#include "model.h"
void init_diff_const(void) {
} /* init_diff_const */

double reac(double* S, int l) {
 double A;
 double B;
 double D;
  A=S[0];
  B=S[1];
  D=S[2];
  switch (l) {
  case 0: return ((-1.000000e+00*(((3.000000e+00*A)*(B*B))-(4.000000e+00
*B))));
break;
  case 1: return (((-1.000000e+00*(((3.000000e+00*A)*(B*B))-(4.000000e+00
*B)))+(-1.000000e+00*(2.000000e+00*B))));
break;
  case 2: return ((2.000000e+00*B));
break;
  } /* switch */
} /* reac */
void eval(double S[equa], double F[equa]) {
 double A;
 double B;
 double D;
  if (S[0]<0.0) S[0]=0.0;
  A=S[0];
  if (S[1]<0.0) S[1]=0.0;
  B=S[1];
  if (S[2]<0.0) S[2]=0.0;
  D=S[2];
  F[0] = (-1.000000e+00*(((3.000000e+00*A)*(B*B))-(4.000000e+00
*B)));
  F[1] = ((-1.000000e+00*(((3.000000e+00*A)*(B*B))-(4.000000e+00
*B)))+(-1.000000e+00*(2.000000e+00*B)));
  F[2] = (2.000000e+00*B);
} /* eval */
void calc_jac(double S[equa], double tk, double R[equa], double Jac[equa+1][equa+1]) 
{
 double A;
 double B;
 double D;
  if (S[0]<0.0) S[0]=0.0;
  A=S[0];
  if (S[1]<0.0) S[1]=0.0;
  B=S[1];
  if (S[2]<0.0) S[2]=0.0;
  D=S[2];
Jac[1][1]=tk*((-1.000000e+00*(3.000000e+00*(B*B))))-(1+4*R[0]);
Jac[1][2]=tk*((-1.000000e+00*(((3.000000e+00*A)*(B*2.000000e+00
))-4.000000e+00)));
Jac[1][3]=tk*(0.000000e+00);
Jac[2][1]=tk*((-1.000000e+00*(3.000000e+00*(B*B))));
Jac[2][2]=tk*(((-1.000000e+00*(((3.000000e+00*A)*(B*2.000000e+00
))-4.000000e+00))+-2.000000e+00))-(1+4*R[1]);
Jac[2][3]=tk*(0.000000e+00);
Jac[3][1]=tk*(0.000000e+00);
Jac[3][2]=tk*(2.000000e+00);
Jac[3][3]=tk*(0.000000e+00)-(1+4*R[2]);
}
reac_F(double S[equa], double K[equa], double R[equa], double tk, int l) {
 double A;
 double B;
 double D;
  A=S[0];
  B=S[1];
  D=S[2];
switch (l) {
case 1:
return (tk*((-1.000000e+00*(((3.000000e+00*A)*(B*B))-(4.000000e+00
*B))))-S[0]*(1+4*R[0])+K[0]);
break;
case 2:
return (tk*(((-1.000000e+00*(((3.000000e+00*A)*(B*B))-(4.000000e+00
*B)))+(-1.000000e+00*(2.000000e+00*B))))-S[1]*(1+4*R[1])+K[1]);
break;
case 3:
return (tk*((2.000000e+00*B))-S[2]*(1+4*R[2])+K[2]);
break;
}
}
void calc_jac2(double S[equa], double Jac[equa][equa]) 
{
double A;
double B;
double D;
A=S[0];
B=S[1];
D=S[2];
Jac[0][0]=(-1.000000e+00*(3.000000e+00*(B*B)));
Jac[0][1]=(-1.000000e+00*(((3.000000e+00*A)*(B*2.000000e+00
))-4.000000e+00));
Jac[0][2]=0.000000e+00;
Jac[1][0]=(-1.000000e+00*(3.000000e+00*(B*B)));
Jac[1][1]=((-1.000000e+00*(((3.000000e+00*A)*(B*2.000000e+00
))-4.000000e+00))+-2.000000e+00);
Jac[1][2]=0.000000e+00;
Jac[2][0]=0.000000e+00;
Jac[2][1]=2.000000e+00;
Jac[2][2]=0.000000e+00;
}
