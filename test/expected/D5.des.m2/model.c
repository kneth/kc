/* ########################
*/

#include <math.h>
#include "model.h"
void init_diff_const(void) {
} /* init_diff_const */

double reac(double* S, int l) {
 double y1;
 double y2;
  y1=S[1];
  y2=S[2];
  switch (l) {
  case 0: return ((1.000000e-02-((1.000000e+00+((y1+1.000000e+03
)*(y1+1.000000e+00)))*((1.000000e-02+y1)+y2))));
break;
  case 1: return ((1.000000e-02-((1.000000e+00+(y2*y2))*((1.000000e-02
+y1)+y2))));
break;
  } /* switch */
} /* reac */
void eval(double S[equa], double F[equa]) {
 double y1;
 double y2;
  y1=S[0];
  y2=S[1];
  F[0] = (1.000000e-02-((1.000000e+00+((y1+1.000000e+03
)*(y1+1.000000e+00)))*((1.000000e-02+y1)+y2)));
  F[1] = (1.000000e-02-((1.000000e+00+(y2*y2))*((1.000000e-02
+y1)+y2)));
} /* eval */
void calc_jac(double S[equa], double tk, double R[equa], double Jac[equa+1][equa+1]) 
{
 double y1;
 double y2;
  y1=S[0];
  y2=S[1];
Jac[1][1]=tk*((-(((((y1+1.000000e+00)+(y1+1.000000e+03))*((
1.000000e-02+y1)+y2))+(1.000000e+00+((y1+1.000000e+03
)*(y1+1.000000e+00)))))))-(1+4*R[0]);
Jac[1][2]=tk*((-((1.000000e+00+((y1+1.000000e+03)*(y1+1.000000e+00
))))));
Jac[2][1]=tk*((-((1.000000e+00+(y2*y2)))));
Jac[2][2]=tk*((-((((y2*2.000000e+00)*((1.000000e-02+y1)+y2)
)+(1.000000e+00+(y2*y2))))))-(1+4*R[1]);
}
double reac_F(double S[equa], double K[equa], double R[equa], double tk, int l) {
 double y1;
 double y2;
  y1=S[0];
  y2=S[1];
switch (l) {
case 1:
return (tk*((1.000000e-02-((1.000000e+00+((y1+1.000000e+03
)*(y1+1.000000e+00)))*((1.000000e-02+y1)+y2))))-S[0]*(1+4*R[0])+K[0]);
break;
case 2:
return (tk*((1.000000e-02-((1.000000e+00+(y2*y2))*((1.000000e-02
+y1)+y2))))-S[1]*(1+4*R[1])+K[1]);
break;
}
}
void calc_jac2(double S[equa], double Jac[equa][equa]) 
{
double y1;
double y2;
y1=S[0];
y2=S[1];
Jac[0][0]=(-(((((y1+1.000000e+00)+(y1+1.000000e+03))*((
1.000000e-02+y1)+y2))+(1.000000e+00+((y1+1.000000e+03
)*(y1+1.000000e+00))))));
Jac[0][1]=(-((1.000000e+00+((y1+1.000000e+03)*(y1+1.000000e+00
)))));
Jac[1][0]=(-((1.000000e+00+(y2*y2))));
Jac[1][1]=(-((((y2*2.000000e+00)*((1.000000e-02+y1)+y2)
)+(1.000000e+00+(y2*y2)))));
}
