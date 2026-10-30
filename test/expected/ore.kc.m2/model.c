/* ########################
*/

#include <math.h>
#include "model.h"
void init_diff_const(void) {
} /* init_diff_const */

double reac(double* S, int l) {
 double HBrO2;
 double P;
 double Br_m;
 double Ce_pppp;
 double BrO3_m;
 double H;
 double HOBr;
 double BrMA;
  HBrO2=S[0];
  Br_m=S[1];
  Ce_pppp=S[2];
  P=0.000000e+00;
  BrO3_m=1.200000e-02;
  H=1.050000e+00;
  HOBr=1.000000e+00;
  BrMA=1.000000e+00;
  switch (l) {
  case 0: return ((((((-1.000000e+00*(1.105000e+00*HBrO2))+((2.400000e-02
*Br_m)*1.102500e+00))+(-1.000000e+00*(((3.000000e+06
*HBrO2)*Br_m)*1.050000e+00)))+(-1.000000e+00*((5.040000e-01
*HBrO2)*1.050000e+00)))+(-2.000000e+00*(3.000000e+03
*(HBrO2*HBrO2)))));
break;
  case 1: return (((((-1.000000e+00*(1.105000e+00*Br_m))+(-1.000000e+00
*((2.400000e-02*Br_m)*1.102500e+00)))+(-1.000000e+00
*(((3.000000e+06*HBrO2)*Br_m)*1.050000e+00)))+(1.320000e-01
*Ce_pppp)));
break;
  case 2: return (((((-1.000000e+00*(1.105000e+00*Ce_pppp))+(2.000000e+00
*((5.040000e-01*HBrO2)*1.050000e+00)))+(-1.000000e+00
*(1.320000e-01*Ce_pppp)))+(-1.000000e+00*(3.500000e-02
*Ce_pppp))));
break;
  } /* switch */
} /* reac */
void eval(double S[equa], double F[equa]) {
 double HBrO2;
 double P;
 double Br_m;
 double Ce_pppp;
 double BrO3_m;
 double H;
 double HOBr;
 double BrMA;
  if (S[0]<0.0) S[0]=0.0;
  HBrO2=S[0];
  if (S[1]<0.0) S[1]=0.0;
  Br_m=S[1];
  if (S[2]<0.0) S[2]=0.0;
  Ce_pppp=S[2];
  P=0.000000e+00;
  BrO3_m=1.200000e-02;
  H=1.050000e+00;
  HOBr=1.000000e+00;
  BrMA=1.000000e+00;
  F[0] = (((((-1.000000e+00*(1.105000e+00*HBrO2))+((2.400000e-02
*Br_m)*1.102500e+00))+(-1.000000e+00*(((3.000000e+06
*HBrO2)*Br_m)*1.050000e+00)))+(-1.000000e+00*((5.040000e-01
*HBrO2)*1.050000e+00)))+(-2.000000e+00*(3.000000e+03
*(HBrO2*HBrO2))));
  F[1] = ((((-1.000000e+00*(1.105000e+00*Br_m))+(-1.000000e+00
*((2.400000e-02*Br_m)*1.102500e+00)))+(-1.000000e+00
*(((3.000000e+06*HBrO2)*Br_m)*1.050000e+00)))+(1.320000e-01
*Ce_pppp));
  F[2] = ((((-1.000000e+00*(1.105000e+00*Ce_pppp))+(2.000000e+00
*((5.040000e-01*HBrO2)*1.050000e+00)))+(-1.000000e+00
*(1.320000e-01*Ce_pppp)))+(-1.000000e+00*(3.500000e-02
*Ce_pppp)));
} /* eval */
void calc_jac(double S[equa], double tk, double R[equa], double Jac[equa+1][equa+1]) 
{
 double HBrO2;
 double P;
 double Br_m;
 double Ce_pppp;
 double BrO3_m;
 double H;
 double HOBr;
 double BrMA;
  if (S[0]<0.0) S[0]=0.0;
  HBrO2=S[0];
  if (S[1]<0.0) S[1]=0.0;
  Br_m=S[1];
  if (S[2]<0.0) S[2]=0.0;
  Ce_pppp=S[2];
  P=0.000000e+00;
  BrO3_m=1.200000e-02;
  H=1.050000e+00;
  HOBr=1.000000e+00;
  BrMA=1.000000e+00;
Jac[1][1]=tk*((((-1.105000e+00+(-1.000000e+00*((3.000000e+06
*Br_m)*1.050000e+00)))+-5.292000e-01)+(-2.000000e+00*
(3.000000e+03*(HBrO2*2.000000e+00)))))-(1+4*R[0]);
Jac[1][2]=tk*((2.646000e-02+(-1.000000e+00*((3.000000e+06
*HBrO2)*1.050000e+00))));
Jac[1][3]=tk*(0.000000e+00);
Jac[2][1]=tk*((-1.000000e+00*((3.000000e+06*Br_m)*1.050000e+00
)));
Jac[2][2]=tk*((-1.131460e+00+(-1.000000e+00*((3.000000e+06
*HBrO2)*1.050000e+00))))-(1+4*R[1]);
Jac[2][3]=tk*(1.320000e-01);
Jac[3][1]=tk*(1.058400e+00);
Jac[3][2]=tk*(0.000000e+00);
Jac[3][3]=tk*(-1.272000e+00)-(1+4*R[2]);
}
double reac_F(double S[equa], double K[equa], double R[equa], double tk, int l) {
 double HBrO2;
 double P;
 double Br_m;
 double Ce_pppp;
 double BrO3_m;
 double H;
 double HOBr;
 double BrMA;
  HBrO2=S[0];
  Br_m=S[1];
  Ce_pppp=S[2];
  P=0.000000e+00;
  BrO3_m=1.200000e-02;
  H=1.050000e+00;
  HOBr=1.000000e+00;
  BrMA=1.000000e+00;
switch (l) {
case 1:
return (tk*((((((-1.000000e+00*(1.105000e+00*HBrO2))+((2.400000e-02
*Br_m)*1.102500e+00))+(-1.000000e+00*(((3.000000e+06
*HBrO2)*Br_m)*1.050000e+00)))+(-1.000000e+00*((5.040000e-01
*HBrO2)*1.050000e+00)))+(-2.000000e+00*(3.000000e+03
*(HBrO2*HBrO2)))))-S[0]*(1+4*R[0])+K[0]);
break;
case 2:
return (tk*(((((-1.000000e+00*(1.105000e+00*Br_m))+(-1.000000e+00
*((2.400000e-02*Br_m)*1.102500e+00)))+(-1.000000e+00
*(((3.000000e+06*HBrO2)*Br_m)*1.050000e+00)))+(1.320000e-01
*Ce_pppp)))-S[1]*(1+4*R[1])+K[1]);
break;
case 3:
return (tk*(((((-1.000000e+00*(1.105000e+00*Ce_pppp))+(2.000000e+00
*((5.040000e-01*HBrO2)*1.050000e+00)))+(-1.000000e+00
*(1.320000e-01*Ce_pppp)))+(-1.000000e+00*(3.500000e-02
*Ce_pppp))))-S[2]*(1+4*R[2])+K[2]);
break;
}
}
void calc_jac2(double S[equa], double Jac[equa][equa]) 
{
double HBrO2;
double P;
double Br_m;
double Ce_pppp;
double BrO3_m;
double H;
double HOBr;
double BrMA;
HBrO2=S[0];
Br_m=S[1];
Ce_pppp=S[2];
  P=0.000000e+00;
  BrO3_m=1.200000e-02;
  H=1.050000e+00;
  HOBr=1.000000e+00;
  BrMA=1.000000e+00;
Jac[0][0]=(((-1.105000e+00+(-1.000000e+00*((3.000000e+06
*Br_m)*1.050000e+00)))+-5.292000e-01)+(-2.000000e+00*
(3.000000e+03*(HBrO2*2.000000e+00))));
Jac[0][1]=(2.646000e-02+(-1.000000e+00*((3.000000e+06
*HBrO2)*1.050000e+00)));
Jac[0][2]=0.000000e+00;
Jac[1][0]=(-1.000000e+00*((3.000000e+06*Br_m)*1.050000e+00
));
Jac[1][1]=(-1.131460e+00+(-1.000000e+00*((3.000000e+06
*HBrO2)*1.050000e+00)));
Jac[1][2]=1.320000e-01;
Jac[2][0]=1.058400e+00;
Jac[2][1]=0.000000e+00;
Jac[2][2]=-1.272000e+00;
}
