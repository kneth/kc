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
 double y5;
 double y6;
  y1=S[1];
  y2=S[2];
  y3=S[3];
  y4=S[4];
  y5=S[5];
  y6=S[6];
  switch (l) {
  case 0: return (((-1.000000e+01*y1)+(1.000000e+02*y2)));
break;
  case 1: return (((-1.000000e+02*y1)-(1.000000e+01*y2)));
break;
  case 2: return ((-4.000000e+00*y3));
break;
  case 3: return ((-(y4)));
break;
  case 4: return ((-5.000000e-01*y5));
break;
  case 5: return ((-1.000000e-01*y6));
break;
  } /* switch */
} /* reac */
void eval(double S[equa], double F[equa]) {
 double y1;
 double y2;
 double y3;
 double y4;
 double y5;
 double y6;
  y1=S[0];
  y2=S[1];
  y3=S[2];
  y4=S[3];
  y5=S[4];
  y6=S[5];
  F[0] = ((-1.000000e+01*y1)+(1.000000e+02*y2));
  F[1] = ((-1.000000e+02*y1)-(1.000000e+01*y2));
  F[2] = (-4.000000e+00*y3);
  F[3] = (-(y4));
  F[4] = (-5.000000e-01*y5);
  F[5] = (-1.000000e-01*y6);
} /* eval */
void calc_jac(double S[equa], double tk, double R[equa], double Jac[equa+1][equa+1]) 
{
 double y1;
 double y2;
 double y3;
 double y4;
 double y5;
 double y6;
  y1=S[0];
  y2=S[1];
  y3=S[2];
  y4=S[3];
  y5=S[4];
  y6=S[5];
Jac[1][1]=tk*(-1.000000e+01)-(1+4*R[0]);
Jac[1][2]=tk*(1.000000e+02);
Jac[1][3]=tk*(0.000000e+00);
Jac[1][4]=tk*(0.000000e+00);
Jac[1][5]=tk*(0.000000e+00);
Jac[1][6]=tk*(0.000000e+00);
Jac[2][1]=tk*(-1.000000e+02);
Jac[2][2]=tk*(-1.000000e+01)-(1+4*R[1]);
Jac[2][3]=tk*(0.000000e+00);
Jac[2][4]=tk*(0.000000e+00);
Jac[2][5]=tk*(0.000000e+00);
Jac[2][6]=tk*(0.000000e+00);
Jac[3][1]=tk*(0.000000e+00);
Jac[3][2]=tk*(0.000000e+00);
Jac[3][3]=tk*(-4.000000e+00)-(1+4*R[2]);
Jac[3][4]=tk*(0.000000e+00);
Jac[3][5]=tk*(0.000000e+00);
Jac[3][6]=tk*(0.000000e+00);
Jac[4][1]=tk*(-0.000000e+00);
Jac[4][2]=tk*(-0.000000e+00);
Jac[4][3]=tk*(-0.000000e+00);
Jac[4][4]=tk*(-1.000000e+00)-(1+4*R[3]);
Jac[4][5]=tk*(-0.000000e+00);
Jac[4][6]=tk*(-0.000000e+00);
Jac[5][1]=tk*(0.000000e+00);
Jac[5][2]=tk*(0.000000e+00);
Jac[5][3]=tk*(0.000000e+00);
Jac[5][4]=tk*(0.000000e+00);
Jac[5][5]=tk*(-5.000000e-01)-(1+4*R[4]);
Jac[5][6]=tk*(0.000000e+00);
Jac[6][1]=tk*(0.000000e+00);
Jac[6][2]=tk*(0.000000e+00);
Jac[6][3]=tk*(0.000000e+00);
Jac[6][4]=tk*(0.000000e+00);
Jac[6][5]=tk*(0.000000e+00);
Jac[6][6]=tk*(-1.000000e-01)-(1+4*R[5]);
}
double reac_F(double S[equa], double K[equa], double R[equa], double tk, int l) {
 double y1;
 double y2;
 double y3;
 double y4;
 double y5;
 double y6;
  y1=S[0];
  y2=S[1];
  y3=S[2];
  y4=S[3];
  y5=S[4];
  y6=S[5];
switch (l) {
case 1:
return (tk*(((-1.000000e+01*y1)+(1.000000e+02*y2)))-S[0]*(1+4*R[0])+K[0]);
break;
case 2:
return (tk*(((-1.000000e+02*y1)-(1.000000e+01*y2)))-S[1]*(1+4*R[1])+K[1]);
break;
case 3:
return (tk*((-4.000000e+00*y3))-S[2]*(1+4*R[2])+K[2]);
break;
case 4:
return (tk*((-(y4)))-S[3]*(1+4*R[3])+K[3]);
break;
case 5:
return (tk*((-5.000000e-01*y5))-S[4]*(1+4*R[4])+K[4]);
break;
case 6:
return (tk*((-1.000000e-01*y6))-S[5]*(1+4*R[5])+K[5]);
break;
}
}
void calc_jac2(double S[equa], double Jac[equa][equa]) 
{
double y1;
double y2;
double y3;
double y4;
double y5;
double y6;
y1=S[0];
y2=S[1];
y3=S[2];
y4=S[3];
y5=S[4];
y6=S[5];
Jac[0][0]=-1.000000e+01;
Jac[0][1]=1.000000e+02;
Jac[0][2]=0.000000e+00;
Jac[0][3]=0.000000e+00;
Jac[0][4]=0.000000e+00;
Jac[0][5]=0.000000e+00;
Jac[1][0]=-1.000000e+02;
Jac[1][1]=-1.000000e+01;
Jac[1][2]=0.000000e+00;
Jac[1][3]=0.000000e+00;
Jac[1][4]=0.000000e+00;
Jac[1][5]=0.000000e+00;
Jac[2][0]=0.000000e+00;
Jac[2][1]=0.000000e+00;
Jac[2][2]=-4.000000e+00;
Jac[2][3]=0.000000e+00;
Jac[2][4]=0.000000e+00;
Jac[2][5]=0.000000e+00;
Jac[3][0]=-0.000000e+00;
Jac[3][1]=-0.000000e+00;
Jac[3][2]=-0.000000e+00;
Jac[3][3]=-1.000000e+00;
Jac[3][4]=-0.000000e+00;
Jac[3][5]=-0.000000e+00;
Jac[4][0]=0.000000e+00;
Jac[4][1]=0.000000e+00;
Jac[4][2]=0.000000e+00;
Jac[4][3]=0.000000e+00;
Jac[4][4]=-5.000000e-01;
Jac[4][5]=0.000000e+00;
Jac[5][0]=0.000000e+00;
Jac[5][1]=0.000000e+00;
Jac[5][2]=0.000000e+00;
Jac[5][3]=0.000000e+00;
Jac[5][4]=0.000000e+00;
Jac[5][5]=-1.000000e-01;
}
