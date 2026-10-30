/* ########################
*/

#include <math.h>
#include "model.h"
void init_diff_const(void) {
} /* init_diff_const */

double reac(double* S, int l) {
 double BrO3_m;
 double Br_m;
 double H_p;
 double HBrO2;
 double HOBr;
 double Ce_pppp;
 double Ce_ppp;
 double P;
  Br_m=S[0];
  HBrO2=S[1];
  Ce_pppp=S[2];
  BrO3_m=1.201089e-02;
  H_p=1.100000e+00;
  HOBr=0.000000e+00;
  Ce_ppp=0.000000e+00;
  P=0.000000e+00;
  switch (l) {
  case 0: return ((((-1.000000e+00*((2.402178e-02*Br_m)*1.210000e+00
))+(-1.000000e+00*(((3.000000e+06*HBrO2)*Br_m)*1.100000e+00
)))+(1.320000e-01*Ce_pppp)));
break;
  case 1: return ((((((2.402178e-02*Br_m)*1.210000e+00)+(-1.000000e+00
*(((3.000000e+06*HBrO2)*Br_m)*1.100000e+00)))+(-2.000000e+00
*(3.000000e+03*(HBrO2*HBrO2))))+((5.044574e-01*HBrO2)*
1.100000e+00)));
break;
  case 2: return ((((2.000000e+00*((5.044574e-01*HBrO2)*1.100000e+00
))+(-1.000000e+00*(1.320000e-01*Ce_pppp)))+(-1.000000e+00
*(3.499000e-02*Ce_pppp))));
break;
  } /* switch */
} /* reac */
void eval(double S[equa], double F[equa]) {
 double BrO3_m;
 double Br_m;
 double H_p;
 double HBrO2;
 double HOBr;
 double Ce_pppp;
 double Ce_ppp;
 double P;
  if (S[0]<0.0) S[0]=0.0;
  Br_m=S[0];
  if (S[1]<0.0) S[1]=0.0;
  HBrO2=S[1];
  if (S[2]<0.0) S[2]=0.0;
  Ce_pppp=S[2];
  BrO3_m=1.201089e-02;
  H_p=1.100000e+00;
  HOBr=0.000000e+00;
  Ce_ppp=0.000000e+00;
  P=0.000000e+00;
  F[0] = (((-1.000000e+00*((2.402178e-02*Br_m)*1.210000e+00
))+(-1.000000e+00*(((3.000000e+06*HBrO2)*Br_m)*1.100000e+00
)))+(1.320000e-01*Ce_pppp));
  F[1] = (((((2.402178e-02*Br_m)*1.210000e+00)+(-1.000000e+00
*(((3.000000e+06*HBrO2)*Br_m)*1.100000e+00)))+(-2.000000e+00
*(3.000000e+03*(HBrO2*HBrO2))))+((5.044574e-01*HBrO2)*
1.100000e+00));
  F[2] = (((2.000000e+00*((5.044574e-01*HBrO2)*1.100000e+00
))+(-1.000000e+00*(1.320000e-01*Ce_pppp)))+(-1.000000e+00
*(3.499000e-02*Ce_pppp)));
} /* eval */
void calc_jac(double S[equa], double tk, double R[equa], double Jac[equa+1][equa+1]) 
{
 double BrO3_m;
 double Br_m;
 double H_p;
 double HBrO2;
 double HOBr;
 double Ce_pppp;
 double Ce_ppp;
 double P;
  if (S[0]<0.0) S[0]=0.0;
  Br_m=S[0];
  if (S[1]<0.0) S[1]=0.0;
  HBrO2=S[1];
  if (S[2]<0.0) S[2]=0.0;
  Ce_pppp=S[2];
  BrO3_m=1.201089e-02;
  H_p=1.100000e+00;
  HOBr=0.000000e+00;
  Ce_ppp=0.000000e+00;
  P=0.000000e+00;
Jac[1][1]=tk*((-2.906636e-02+(-1.000000e+00*((3.000000e+06
*HBrO2)*1.100000e+00))))-(1+4*R[0]);
Jac[1][2]=tk*((-1.000000e+00*((3.000000e+06*Br_m)*1.100000e+00
)));
Jac[1][3]=tk*(1.320000e-01);
Jac[2][1]=tk*((2.906636e-02+(-1.000000e+00*((3.000000e+06
*HBrO2)*1.100000e+00))));
Jac[2][2]=tk*((((-1.000000e+00*((3.000000e+06*Br_m)*1.100000e+00
))+(-2.000000e+00*(3.000000e+03*(HBrO2*2.000000e+00)
)))+5.549031e-01))-(1+4*R[1]);
Jac[2][3]=tk*(0.000000e+00);
Jac[3][1]=tk*(0.000000e+00);
Jac[3][2]=tk*(1.109806e+00);
Jac[3][3]=tk*(-1.669900e-01)-(1+4*R[2]);
}
reac_F(double S[equa], double K[equa], double R[equa], double tk, int l) {
 double BrO3_m;
 double Br_m;
 double H_p;
 double HBrO2;
 double HOBr;
 double Ce_pppp;
 double Ce_ppp;
 double P;
  Br_m=S[0];
  HBrO2=S[1];
  Ce_pppp=S[2];
  BrO3_m=1.201089e-02;
  H_p=1.100000e+00;
  HOBr=0.000000e+00;
  Ce_ppp=0.000000e+00;
  P=0.000000e+00;
switch (l) {
case 1:
return (tk*((((-1.000000e+00*((2.402178e-02*Br_m)*1.210000e+00
))+(-1.000000e+00*(((3.000000e+06*HBrO2)*Br_m)*1.100000e+00
)))+(1.320000e-01*Ce_pppp)))-S[0]*(1+4*R[0])+K[0]);
break;
case 2:
return (tk*((((((2.402178e-02*Br_m)*1.210000e+00)+(-1.000000e+00
*(((3.000000e+06*HBrO2)*Br_m)*1.100000e+00)))+(-2.000000e+00
*(3.000000e+03*(HBrO2*HBrO2))))+((5.044574e-01*HBrO2)*
1.100000e+00)))-S[1]*(1+4*R[1])+K[1]);
break;
case 3:
return (tk*((((2.000000e+00*((5.044574e-01*HBrO2)*1.100000e+00
))+(-1.000000e+00*(1.320000e-01*Ce_pppp)))+(-1.000000e+00
*(3.499000e-02*Ce_pppp))))-S[2]*(1+4*R[2])+K[2]);
break;
}
}
void calc_jac2(double S[equa], double Jac[equa][equa]) 
{
double BrO3_m;
double Br_m;
double H_p;
double HBrO2;
double HOBr;
double Ce_pppp;
double Ce_ppp;
double P;
Br_m=S[0];
HBrO2=S[1];
Ce_pppp=S[2];
  BrO3_m=1.201089e-02;
  H_p=1.100000e+00;
  HOBr=0.000000e+00;
  Ce_ppp=0.000000e+00;
  P=0.000000e+00;
Jac[0][0]=(-2.906636e-02+(-1.000000e+00*((3.000000e+06
*HBrO2)*1.100000e+00)));
Jac[0][1]=(-1.000000e+00*((3.000000e+06*Br_m)*1.100000e+00
));
Jac[0][2]=1.320000e-01;
Jac[1][0]=(2.906636e-02+(-1.000000e+00*((3.000000e+06
*HBrO2)*1.100000e+00)));
Jac[1][1]=(((-1.000000e+00*((3.000000e+06*Br_m)*1.100000e+00
))+(-2.000000e+00*(3.000000e+03*(HBrO2*2.000000e+00)
)))+5.549031e-01);
Jac[1][2]=0.000000e+00;
Jac[2][0]=0.000000e+00;
Jac[2][1]=1.109806e+00;
Jac[2][2]=-1.669900e-01;
}
