#include <stdio.h>
#include <math.h>
#include <conio.h>


double Argu(double g, double h) {
/* return (360.0/2.0/3.14159265359*atan(h,g)) */
   return (57.2957795131*atan2(h,g));
} /* Argu */

double Radius(double a, double b) {
   return (sqrt(a*a+b*b));
} /* Radius */

void compamppha(int n, double P[4][4], double amp[4], double phase[4]) {

   int   i;

   for(i=0;i<n;i++) {
     amp[i]= Radius(P[i][0],P[i][1]);
     phase[i]= Argu(P[i][0],P[i][1]);
   }
} /* compamppha */

void stopdata(int n, int m, double IP[4][4], double x[4],
	      double q[4], double fi[4], double *qd, double *fid) {

   int     i;
   double  gd, hd;

   for(i=0;i<n;i++) {
     q[i]=  1.0/Radius(IP[0][i],IP[1][i]);
     fi[i]= Argu(-IP[0][i],-IP[1][i]);
   }

   gd= 0.0; hd= 0.0;
   for(i=0;i<n;i++) {
     gd += (IP[0][i]*x[i]);
     hd += (IP[1][i]*x[i]);
   }

   gd *= -1.0; hd *= -1.0;
   *qd=  -x[m-1]/Radius(gd,hd);
   *fid= Argu(-gd,-hd);
} /* stopdata */


int main(void)
{

   int    i;
   double ph[4],am[4],x[4],qu[4],fu[4],qd,fd;
   double IP[4][4],P[4][4];

   clrscr();

   P[0][0]=  0.1474;P[0][1]= -0.06801;P[0][2]= -0.04463;P[0][3]= -0.34473;
   P[1][0]= -0.2052;P[1][1]=  0.27627;P[1][2]= -0.04498;P[1][3]= -1.89213;
   P[2][0]=  1.0000;P[2][1]=  0.00000;P[2][2]=  1.00000;P[2][3]=  1.00000;
   P[3][0]=  0.4706;P[3][1]= -0.00665;P[3][2]= -0.53426;P[3][3]=  4.31643;


   IP[0][0]=  3.2273;IP[0][1]=  0.80668;IP[0][2]=  0.45111;IP[0][3]=  0.50685;
   IP[1][0]= -2.6220;IP[1][1]=  2.99771;IP[1][2]=  0.54101;IP[1][3]=  0.97930;
   IP[2][0]= -2.5551;IP[2][1]= -0.64367;IP[2][2]=  0.53146;IP[2][3]= -0.60935;
   IP[3][0]= -0.6722;IP[3][1]= -0.16301;IP[3][2]=  0.01743;IP[3][3]=  0.10250;


   x[0]= 0.02810E-6;x[1]= 0.20968E-6;x[2]= 0.19058E-6;x[3]= 0.12369E-6;

   compamppha(4, P, am, ph);
   for(i=0;i<4;i++)
      printf("Amp: %lf .  Phase: %lf\n", am[i], ph[i]);
   printf("\n");

   stopdata(4,3,IP,x,qu,fu,&qd,&fd);
   for(i=0;i<4;i++)
      printf("Que: %lf .  Phase: %lf\n", qu[i], fu[i]);
   printf("Dqu: %lf . Dphas: %lf\n\n",qd,fd);

   return(0);
}

