/* ########################
*/

#include <math.h>
#include "model.h"
void init_diff_const(void) {
} /* init_diff_const */

double reac(double* S, int l) {
 double u1;
 double u2;
 double u3;
 double u4;
  u1=S[1];
  u2=S[2];
  u3=S[3];
  u4=S[4];
  switch (l) {
  case 0: return (log(u1));
break;
  case 1: return (log10(u2));
break;
  case 2: return (((exp(u3)*sin(u3))*cos(u3)));
break;
  case 3: return ((pow(u1,u2)));
break;
  } /* switch */
} /* reac */
void eval(double S[equa], double F[equa]) {
 double u1;
 double u2;
 double u3;
 double u4;
  u1=S[0];
  u2=S[1];
  u3=S[2];
  u4=S[3];
  F[0] = log(u1);
  F[1] = log10(u2);
  F[2] = ((exp(u3)*sin(u3))*cos(u3));
  F[3] = (pow(u1,u2));
} /* eval */
void calc_jac(double S[equa], double tk, double R[equa], double Jac[equa+1][equa+1]) 
{
 double u1;
 double u2;
 double u3;
 double u4;
  u1=S[0];
  u2=S[1];
  u3=S[2];
  u4=S[3];
Jac[1][1]=tk*((1.000000e+00/u1))-(1+4*R[0]);
Jac[1][2]=tk*(0.000000e+00);
Jac[1][3]=tk*(0.000000e+00);
Jac[1][4]=tk*(0.000000e+00);
Jac[2][1]=tk*(0.000000e+00);
Jac[2][2]=tk*((4.342945e-01/u2))-(1+4*R[1]);
Jac[2][3]=tk*(0.000000e+00);
Jac[2][4]=tk*(0.000000e+00);
Jac[3][1]=tk*(0.000000e+00);
Jac[3][2]=tk*(0.000000e+00);
Jac[3][3]=tk*(((((exp(u3)*sin(u3))+(exp(u3)*cos(u3)))*cos(u3))+(
(exp(u3)*sin(u3))*(-(sin(u3))))))-(1+4*R[2]);
Jac[3][4]=tk*(0.000000e+00);
Jac[4][1]=tk*(((pow(u1,(u2-1.000000e+00)))*u2));
Jac[4][2]=tk*(((pow(u1,u2))*log(u1)));
Jac[4][3]=tk*(0.000000e+00);
Jac[4][4]=tk*(0.000000e+00)-(1+4*R[3]);
}
double reac_F(double S[equa], double K[equa], double R[equa], double tk, int l) {
 double u1;
 double u2;
 double u3;
 double u4;
  u1=S[0];
  u2=S[1];
  u3=S[2];
  u4=S[3];
switch (l) {
case 1:
return (tk*(log(u1))-S[0]*(1+4*R[0])+K[0]);
break;
case 2:
return (tk*(log10(u2))-S[1]*(1+4*R[1])+K[1]);
break;
case 3:
return (tk*(((exp(u3)*sin(u3))*cos(u3)))-S[2]*(1+4*R[2])+K[2]);
break;
case 4:
return (tk*((pow(u1,u2)))-S[3]*(1+4*R[3])+K[3]);
break;
}
}
void calc_jac2(double S[equa], double Jac[equa][equa]) 
{
double u1;
double u2;
double u3;
double u4;
u1=S[0];
u2=S[1];
u3=S[2];
u4=S[3];
Jac[0][0]=(1.000000e+00/u1);
Jac[0][1]=0.000000e+00;
Jac[0][2]=0.000000e+00;
Jac[0][3]=0.000000e+00;
Jac[1][0]=0.000000e+00;
Jac[1][1]=(4.342945e-01/u2);
Jac[1][2]=0.000000e+00;
Jac[1][3]=0.000000e+00;
Jac[2][0]=0.000000e+00;
Jac[2][1]=0.000000e+00;
Jac[2][2]=((((exp(u3)*sin(u3))+(exp(u3)*cos(u3)))*cos(u3))+(
(exp(u3)*sin(u3))*(-(sin(u3)))));
Jac[2][3]=0.000000e+00;
Jac[3][0]=((pow(u1,(u2-1.000000e+00)))*u2);
Jac[3][1]=((pow(u1,u2))*log(u1));
Jac[3][2]=0.000000e+00;
Jac[3][3]=0.000000e+00;
}
