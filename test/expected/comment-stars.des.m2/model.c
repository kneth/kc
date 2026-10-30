/* ########################
*/

#include <math.h>
#include "model.h"
void init_diff_const(void) {
} /* init_diff_const */

double reac(double* S, int l) {
 double x;
  x=S[1];
  switch (l) {
  case 0: return ((-(x)));
break;
  } /* switch */
} /* reac */
void eval(double S[equa], double F[equa]) {
 double x;
  x=S[0];
  F[0] = (-(x));
} /* eval */
void calc_jac(double S[equa], double tk, double R[equa], double Jac[equa+1][equa+1]) 
{
 double x;
  x=S[0];
Jac[1][1]=tk*(-1.000000e+00)-(1+4*R[0]);
}
double reac_F(double S[equa], double K[equa], double R[equa], double tk, int l) {
 double x;
  x=S[0];
switch (l) {
case 1:
return (tk*((-(x)))-S[0]*(1+4*R[0])+K[0]);
break;
}
}
void calc_jac2(double S[equa], double Jac[equa][equa]) 
{
double x;
x=S[0];
Jac[0][0]=-1.000000e+00;
}
