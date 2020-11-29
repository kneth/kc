#include <stdio.h>
#include <math.h>
#include "quench.h"

#define novar_ 4
#define noreac_ 11
double x_[novar_],v_[novar_];
double jacobi_matx[novar_][novar_];
int ref;
char species[25][novar_];
double rfw_[noreac_];
double rrv_[noreac_];
double rfwds_[noreac_][novar_];
double rrvds_[noreac_][novar_];

void ReacRate(double *S_) {
double HBrO2;
double Brm;
double Ce4p;
double HBrO;
HBrO2=S_[0];
Brm=S_[1];
Ce4p=S_[2];
HBrO=S_[3];
v_[0]=((((((-1.000000e+00*(4.709600e-05*HBrO2))+((2.400000e-02
*Brm)*4.900000e-01))+(-1.000000e+00*(((3.000000e+06*
HBrO2)*Brm)*7.000000e-01)))+((5.040000e-01*HBrO2)*7.000000e-01
))+(-2.000000e+00*(3.000000e+03*(HBrO2*HBrO2)))));

v_[1]=((((((-1.000000e+00*(4.709600e-05*Brm))+(-1.000000e+00
*((2.400000e-02*Brm)*4.900000e-01)))+(-1.000000e+00*
(((3.000000e+06*HBrO2)*Brm)*7.000000e-01)))+(2.500000e-01
*(1.040000e-01*Ce4p)))+(8.000000e-02*HBrO)));

v_[2]=((((-1.000000e+00*(4.709600e-05*Ce4p))+(2.000000e+00
*((5.040000e-01*HBrO2)*7.000000e-01)))+(-1.000000e+00
*(1.040000e-01*Ce4p))));
v_[3]=(((((((-1.000000e+00*(4.709600e-05*HBrO))+((2.400000e-02
*Brm)*4.900000e-01))+(2.000000e+00*(((3.000000e+06*
HBrO2)*Brm)*7.000000e-01)))+(3.000000e+03*(HBrO2*HBrO2
)))+(-1.000000e+00*(8.000000e-02*HBrO)))+(-1.000000e+00
*(1.400000e-01*HBrO))));
}

void jacobi(double *S_){
double HBrO2;
double Brm;
double Ce4p;
double HBrO;
HBrO2=S_[0];
Brm=S_[1];
Ce4p=S_[2];
HBrO=S_[3];
jacobi_matx[0][0]=(((-4.709600e-05+(-1.000000e+00*((3.000000e+06
*Brm)*7.000000e-01)))+3.528000e-01)+(-2.000000e+00*(
3.000000e+03*(HBrO2*2.000000e+00))));
jacobi_matx[0][1]=(1.176000e-02+(-1.000000e+00*((3.000000e+06
*HBrO2)*7.000000e-01)));
jacobi_matx[1][0]=(-1.000000e+00*((3.000000e+06*Brm)*7.000000e-01
));
jacobi_matx[1][1]=(-1.180710e-02+(-1.000000e+00*((3.000000e+06
*HBrO2)*7.000000e-01)));
jacobi_matx[3][0]=((2.000000e+00*((3.000000e+06*Brm)*7.000000e-01
))+(3.000000e+03*(HBrO2*2.000000e+00)));
jacobi_matx[3][1]=(1.176000e-02+(2.000000e+00*((3.000000e+06
*HBrO2)*7.000000e-01)));
}

void InitValues(void) {
x_[0]=2.810000e-08;
x_[1]=2.096800e-07;
x_[2]=1.905800e-07;
x_[3]=1.236900e-07;
(void) strcpy(species[0], "HBrO2");
(void) strcpy(species[1], "Brm");
(void) strcpy(species[2], "Ce4p");
(void) strcpy(species[3], "HBrO");

ref= 3;

jacobi_matx[0][2]=0.000000e+00;
jacobi_matx[0][3]=0.000000e+00;
jacobi_matx[1][2]=2.600000e-02;
jacobi_matx[1][3]=8.000000e-02;
jacobi_matx[2][0]=7.056000e-01;
jacobi_matx[2][1]=0.000000e+00;
jacobi_matx[2][2]=-1.040471e-01;
jacobi_matx[2][3]=0.000000e+00;
jacobi_matx[3][2]=0.000000e+00;
jacobi_matx[3][3]=-2.200471e-01;
}


void ReacFlow(double *S_) {
double HBrO2;
double Brm;
double Ce4p;
double HBrO;
HBrO2=S_[0];
Brm=S_[1];
Ce4p=S_[2];
HBrO=S_[3];

rfw_[0]= 4.709600e-05*HBrO2;                  rrv_[0]= 0.0;
rfw_[1]= 4.709600e-05*Brm;                    rrv_[1]= 0.0;
rfw_[2]= 4.709600e-05*Ce4p;                   rrv_[2]= 0.0;
rfw_[3]= 4.709600e-05*HBrO;                   rrv_[3]= 0.0;

rfw_[4]= 2.400000e-02*Brm*4.900000e-01;       rrv_[4]= 0.0;
rfw_[5]= 3.000000e+06*HBrO2*Brm*7.000000e-01; rrv_[5]= 0.0;
rfw_[6]= 5.040000e-01*HBrO2*7.000000e-01;     rrv_[6]= 0.0;
rfw_[7]= 3.000000e+03*HBrO2*HBrO2;            rrv_[7]= 0.0;
rfw_[8]= 1.040000e-01*Ce4p;                   rrv_[8]= 0.0;
rfw_[9]= 8.000000e-02*HBrO;                   rrv_[9]= 0.0;
rfw_[10]=1.400000e-01*HBrO;                   rrv_[10]=0.0;
}

void ReacFlowdS(double *S_) {
int i,j;
double HBrO2;
double Brm;
double Ce4p;
double HBrO;
HBrO2=S_[0];
Brm=S_[1];
Ce4p=S_[2];
HBrO=S_[3];

for(i=0;i<noreac_;i++)
 for(j=0;j<novar_;j++) {
   rfwds_[i][j]= 0.0;
   rrvds_[i][j]= 0.0;
 }

rfwds_[0][0]= 4.709600e-05;
rfwds_[1][1]= 4.709600e-05;
rfwds_[2][2]= 4.709600e-05;
rfwds_[3][3]= 4.709600e-05;

rfwds_[4][1]= 2.400000e-02*4.900000e-01;
rfwds_[5][0]= 3.000000e+06*Brm*7.000000e-01;
rfwds_[5][1]= 3.000000e+06*HBrO2*7.000000e-01;
rfwds_[6][0]= 5.040000e-01*7.000000e-01;
rfwds_[7][0]= 3.000000e+03*2.0*HBrO2;
rfwds_[8][2]= 1.040000e-01;
rfwds_[9][3]= 8.000000e-02;
rfwds_[10][3]=1.400000e-01;
}

int main(void)
{

   int    i;
   double *ph[novar_],*am[novar_],x[novar_],qu[novar_],fu[novar_],qd,fd;
   double IP[novar_][novar_],P[novar_][novar_];

   InitValues();

   P[0][0]=  0.1474;P[0][1]= -0.06801;P[0][2]= -0.04463;P[0][3]= -0.34473;
   P[1][0]= -0.2052;P[1][1]=  0.27627;P[1][2]= -0.04498;P[1][3]= -1.89213;
   P[2][0]=  1.0000;P[2][1]=  0.00000;P[2][2]=  1.00000;P[2][3]=  1.00000;
   P[3][0]=  0.4706;P[3][1]= -0.00665;P[3][2]= -0.53426;P[3][3]=  4.31643;


   IP[0][0]=  3.2273;IP[0][1]=  0.80668;IP[0][2]=  0.45111;IP[0][3]=  0.50685;
   IP[1][0]= -2.6220;IP[1][1]=  2.99771;IP[1][2]=  0.54101;IP[1][3]=  0.97930;
   IP[2][0]= -2.5551;IP[2][1]= -0.64367;IP[2][2]=  0.53146;IP[2][3]= -0.60935;
   IP[3][0]= -0.6722;IP[3][1]= -0.16301;IP[3][2]=  0.01743;IP[3][3]=  0.10250;


   /* Stationary point */
   ReacRate(x_);
   printf("\nStationary point\n");
   for(i=0;i<novar_;i++)
      printf("%s\t%le\t%le\n",species[i],x_[i],v_[i]);

   /* Eigenvalues and eigenvectors */

   /* Call to hqr2alg, and print out of the eigenvalues and -vectors */
   /* for this example eigenvectors in P */

   /* Inverse matrix of the eigenvector matrix */

   /* Call to InversMat, for this example in IP */

   /* Stationary point: Complex amplitudes and phases */

   printf("\nStationary point\n");
   compamppha(novar_, P, am, ph);
   for(i=0;i<novar_;i++)
      printf("Amp: %10.5lf\tPhase: %6.2lf\n", am[i], ph[i]);
   printf("\n");

   /* Stationary point: Quenching data */

   printf("\nStationary point\n");
   stopdata(novar_,ref,IP,x_,qu,fu,&qd,&fd);
   for(i=0;i<novar_;i++)
      printf("Que: %10.5lf\tPhase: %6.2lf\n", qu[i], fu[i]);
   printf("Dil. Que: %10.5lf\tPhase: %6.2lf\n\n",qd,fd);

   /* Stationary reaction flow */
   printf("Stationary Reaction Flow\n");
   ReacFlow(x_);
   for(i=0;i<noreac_;i++)
      printf("Reac %d: Forward: %10.5le\tReverse: %10.5le\n", i+1,rfw_[i], rrv_[i]);
   printf("\n");



   return(0);
}
