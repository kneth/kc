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
 double y4;
  y1=S[1];
  y2=S[2];
  y3=S[3];
  y4=S[4];
  switch (l) {
  case 0: return ((y3-((1.000000e+02*y1)*y2)));
break;
  case 1: return ((((y3+(2.000000e+00*y4))-((1.000000e+02*y1)*
y2))-(2.000000e+04*(y2*y2))));
break;
  case 2: return (((-(y3))+((1.000000e+02*y1)*y2)));
break;
  case 3: return (((-(y4))+(1.000000e+04*(y2*y2))));
break;
  } /* switch */
} /* reac */
void eval(double S[equa], double F[equa]) {
 double y1;
 double y2;
 double y3;
 double y4;
  y1=S[0];
  y2=S[1];
  y3=S[2];
  y4=S[3];
  F[0] = (y3-((1.000000e+02*y1)*y2));
  F[1] = (((y3+(2.000000e+00*y4))-((1.000000e+02*y1)*
y2))-(2.000000e+04*(y2*y2)));
  F[2] = ((-(y3))+((1.000000e+02*y1)*y2));
  F[3] = ((-(y4))+(1.000000e+04*(y2*y2)));
} /* eval */
void calc_jac(double S[equa], double tk, double R[equa], double Jac[equa+1][equa+1]) 
{
 double y1;
 double y2;
 double y3;
 double y4;
  y1=S[0];
  y2=S[1];
  y3=S[2];
  y4=S[3];
Jac[1][1]=tk*((-((1.000000e+02*y2))))-(1+4*R[0]);
Jac[1][2]=tk*((-((1.000000e+02*y1))));
Jac[1][3]=tk*(1.000000e+00);
Jac[1][4]=tk*(0.000000e+00);
Jac[2][1]=tk*((-((1.000000e+02*y2))));
Jac[2][2]=tk*((-(((1.000000e+02*y1)+(2.000000e+04*(y2*2.000000e+00
))))))-(1+4*R[1]);
Jac[2][3]=tk*(1.000000e+00);
Jac[2][4]=tk*(2.000000e+00);
Jac[3][1]=tk*((1.000000e+02*y2));
Jac[3][2]=tk*((1.000000e+02*y1));
Jac[3][3]=tk*(-1.000000e+00)-(1+4*R[2]);
Jac[3][4]=tk*(0.000000e+00);
Jac[4][1]=tk*(0.000000e+00);
Jac[4][2]=tk*((1.000000e+04*(y2*2.000000e+00)));
Jac[4][3]=tk*(0.000000e+00);
Jac[4][4]=tk*(-1.000000e+00)-(1+4*R[3]);
}
reac_F(double S[equa], double K[equa], double R[equa], double tk, int l) {
 double y1;
 double y2;
 double y3;
 double y4;
  y1=S[0];
  y2=S[1];
  y3=S[2];
  y4=S[3];
switch (l) {
case 1:
return (tk*((y3-((1.000000e+02*y1)*y2)))-S[0]*(1+4*R[0])+K[0]);
break;
case 2:
return (tk*((((y3+(2.000000e+00*y4))-((1.000000e+02*y1)*
y2))-(2.000000e+04*(y2*y2))))-S[1]*(1+4*R[1])+K[1]);
break;
case 3:
return (tk*(((-(y3))+((1.000000e+02*y1)*y2)))-S[2]*(1+4*R[2])+K[2]);
break;
case 4:
return (tk*(((-(y4))+(1.000000e+04*(y2*y2))))-S[3]*(1+4*R[3])+K[3]);
break;
}
}
void calc_jac2(double S[equa], double Jac[equa][equa]) 
{
double y1;
double y2;
double y3;
double y4;
y1=S[0];
y2=S[1];
y3=S[2];
y4=S[3];
Jac[0][0]=(-((1.000000e+02*y2)));
Jac[0][1]=(-((1.000000e+02*y1)));
Jac[0][2]=1.000000e+00;
Jac[0][3]=0.000000e+00;
Jac[1][0]=(-((1.000000e+02*y2)));
Jac[1][1]=(-(((1.000000e+02*y1)+(2.000000e+04*(y2*2.000000e+00
)))));
Jac[1][2]=1.000000e+00;
Jac[1][3]=2.000000e+00;
Jac[2][0]=(1.000000e+02*y2);
Jac[2][1]=(1.000000e+02*y1);
Jac[2][2]=-1.000000e+00;
Jac[2][3]=0.000000e+00;
Jac[3][0]=0.000000e+00;
Jac[3][1]=(1.000000e+04*(y2*2.000000e+00));
Jac[3][2]=0.000000e+00;
Jac[3][3]=-1.000000e+00;
}
