/* ########################
*/

#include <math.h>
#include "model.h"
void init_diff_const(void) {
} /* init_diff_const */

double reac(double* S, int l) {
 double x;
 double y;
 double z;
 double w;
  x=S[1];
  y=S[2];
  z=S[3];
  w=S[4];
  switch (l) {
  case 0: return (((3.000000e+00*x)+(6.000000e+00*y)));
break;
  case 1: return ((((2.000000e+00*(x*y))-(2.000000e+00*(x*z)))
-3.000000e+00));
break;
  case 2: return (((-((((x+(x*y))+(x*y))+(x/y))))+(2.000000e+00*(x/
y))));
break;
  case 3: return (((((-2.000000e+00*x)+x)-(2.000000e+00*y))-(x+
(2.000000e+00*z))));
break;
  } /* switch */
} /* reac */
void eval(double S[equa], double F[equa]) {
 double x;
 double y;
 double z;
 double w;
  x=S[0];
  y=S[1];
  z=S[2];
  w=S[3];
  F[0] = ((3.000000e+00*x)+(6.000000e+00*y));
  F[1] = (((2.000000e+00*(x*y))-(2.000000e+00*(x*z)))
-3.000000e+00);
  F[2] = ((-((((x+(x*y))+(x*y))+(x/y))))+(2.000000e+00*(x/
y)));
  F[3] = ((((-2.000000e+00*x)+x)-(2.000000e+00*y))-(x+
(2.000000e+00*z)));
} /* eval */
void calc_jac(double S[equa], double tk, double R[equa], double Jac[equa+1][equa+1]) 
{
 double x;
 double y;
 double z;
 double w;
  x=S[0];
  y=S[1];
  z=S[2];
  w=S[3];
Jac[1][1]=tk*(3.000000e+00)-(1+4*R[0]);
Jac[1][2]=tk*(6.000000e+00);
Jac[1][3]=tk*(0.000000e+00);
Jac[1][4]=tk*(0.000000e+00);
Jac[2][1]=tk*(((2.000000e+00*y)-(2.000000e+00*z)));
Jac[2][2]=tk*((2.000000e+00*x))-(1+4*R[1]);
Jac[2][3]=tk*((-2.000000e+00*x));
Jac[2][4]=tk*(0.000000e+00);
Jac[3][1]=tk*(((-((((1.000000e+00+y)+y)+(y/(y*y)))))+(2.000000e+00
*(y/(y*y)))));
Jac[3][2]=tk*((-((((x+x)-(x/(y*y)))+(2.000000e+00*(x/(y*y))))))
);
Jac[3][3]=tk*(0.000000e+00)-(1+4*R[2]);
Jac[3][4]=tk*(0.000000e+00);
Jac[4][1]=tk*(-2.000000e+00);
Jac[4][2]=tk*(-2.000000e+00);
Jac[4][3]=tk*(-2.000000e+00);
Jac[4][4]=tk*(0.000000e+00)-(1+4*R[3]);
}
double reac_F(double S[equa], double K[equa], double R[equa], double tk, int l) {
 double x;
 double y;
 double z;
 double w;
  x=S[0];
  y=S[1];
  z=S[2];
  w=S[3];
switch (l) {
case 1:
return (tk*(((3.000000e+00*x)+(6.000000e+00*y)))-S[0]*(1+4*R[0])+K[0]);
break;
case 2:
return (tk*((((2.000000e+00*(x*y))-(2.000000e+00*(x*z)))
-3.000000e+00))-S[1]*(1+4*R[1])+K[1]);
break;
case 3:
return (tk*(((-((((x+(x*y))+(x*y))+(x/y))))+(2.000000e+00*(x/
y))))-S[2]*(1+4*R[2])+K[2]);
break;
case 4:
return (tk*(((((-2.000000e+00*x)+x)-(2.000000e+00*y))-(x+
(2.000000e+00*z))))-S[3]*(1+4*R[3])+K[3]);
break;
}
}
void calc_jac2(double S[equa], double Jac[equa][equa]) 
{
double x;
double y;
double z;
double w;
x=S[0];
y=S[1];
z=S[2];
w=S[3];
Jac[0][0]=3.000000e+00;
Jac[0][1]=6.000000e+00;
Jac[0][2]=0.000000e+00;
Jac[0][3]=0.000000e+00;
Jac[1][0]=((2.000000e+00*y)-(2.000000e+00*z));
Jac[1][1]=(2.000000e+00*x);
Jac[1][2]=(-2.000000e+00*x);
Jac[1][3]=0.000000e+00;
Jac[2][0]=((-((((1.000000e+00+y)+y)+(y/(y*y)))))+(2.000000e+00
*(y/(y*y))));
Jac[2][1]=(-((((x+x)-(x/(y*y)))+(2.000000e+00*(x/(y*y))))))
;
Jac[2][2]=0.000000e+00;
Jac[2][3]=0.000000e+00;
Jac[3][0]=-2.000000e+00;
Jac[3][1]=-2.000000e+00;
Jac[3][2]=-2.000000e+00;
Jac[3][3]=0.000000e+00;
}
