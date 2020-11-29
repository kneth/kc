/************************************************************************* 
  This program is solving reaction-diffusion in 1 spatial dim. 
  with periodic boundry conditions.

  (C) Copyright 1993-1996
  Kenneth Geisshirt (kneth@osc.kiku.dk)
  Department of Life Sciences and Chemistry
  Roskilde University
  Marbjergvej 35
  4000 Roskilde
  Denmark

  References:
  [1]    Pattern formation in reaction-diffusion systems.
         O. Jensen and V. O. Pannbacker.
  [2]    Numerical Recipes, 2nd ed.
         W. H. Press et al.
  [3]    Matrix Computations.
         G. H. Golub et al.
  [4]    Numerical Methods for Partial Differential Equations
         W. F. Ames.
  [5]    Numerical Analysis
         D. Kincaid et al.
  [6]    Proc. IEEE, April 1968, p. 744
         D. A. Calahan
 
  The program uses a simpel Crank-Nicolson method (see [4], pp.
  50-51) and Gauss-Seidel to solve lin. equations (see [5], pp.
  189-191 and [3], p. 356). A 4th order Runge-Kutta to integrate
  reaction terms (see [2], p. 711) or Calahan's method (see [6]).
  The lin. equations can also be solved by a tricyclic solver (see
  [2], pp. 74-74). 

  Some features are determined by macros, and they follow here:
  o LINSOLVE:
    1  =  Gauss-Seidel,
    2  =  Cyclic.

  Last updated: 26 March 1996
*************************************************************************/

#define LINSOLVE 2

#define EPS       1e-10     /* used in GaussSeidel */
#define MAX_ITER  125       /* used in GaussSeidel */
#define LAMBDA    0.5       /* used by Crack-Nicolson */

#include <stdio.h>
#include <math.h>
#include "model.h"
#include "randgen.h"

double u[n_grids][equa];
double dx, dt, L, max_t;
double R[equa];
char prefix[25];
double ptime1, ptime2, period;


/***************************************************************************
  The following routines are taken from [1], but modified by the
  author.
***************************************************************************/

/* (C) Copr. 1986-92 Numerical Recipes Software +,. */

void tridag(double a[], double b[], double c[], double r[], 
	    double u[], int n) {

  int j;
  double bet, gam[n_grids];
  
  if (b[0] == 0.0) fprintf(stderr, "Error 1 in tridag");
  u[0]=r[0]/(bet=b[0]);
  for (j=2;j<=n_grids;j++) {
    gam[j-1]=c[j-2]/bet;
    bet=b[j-1]-a[j-1]*gam[j-1];
    if (bet == 0.0) fprintf(stderr, "Error 2 in tridag");
    u[j-1]=(r[j-1]-a[j-1]*u[j-2])/bet;
  }
  for (j=(n_grids-1);j>=1;j--)
    u[j-1] -= gam[j]*u[j];
}

void cyclic(double a[], double b[], double c[], 
	    double alpha, double beta, double r[], double x[], 
	    int n) {

  int i;
  double fact, gamma, bb[n_grids], u[n_grids], z[n_grids];
  
  if (n <= 2) fprintf(stderr, "n too small in cyclic\n");
  gamma = -b[0];
  bb[0]=b[0]-gamma;
  bb[n_grids-1]=b[n_grids-1]-alpha*beta/gamma;
  for (i=2;i<n_grids;i++) bb[i-1]=b[i-1];
  tridag(a,bb,c,r,x,n_grids);
  u[0]=gamma;
  u[n_grids-1]=alpha;
  for (i=2;i<n_grids;i++) u[i-1]=0.0;
  tridag(a,bb,c,u,z,n_grids);
  fact=(x[0]+beta*x[n_grids-1]/gamma)/(1.0+z[0]+beta*z[n_grids-1]/gamma);
  for (i=1;i<=n_grids;i++) x[i-1] -= fact*z[i-1];
}

#define TINY 1e-20

void ludcmp(double a[equa][equa], int n, int indx[equa], double *d) {

  int i,imax,j,k;
  double big,dum,sum,temp;
  double vv[equa];

  *d=1.0;
  for (i=1;i<=n;i++) {
    big=0.0;
    for (j=1;j<=n;j++)
      if ((temp=fabs(a[i-1][j-1])) > big) big=temp;
    if (big == 0.0) fprintf(stderr, "Singular matrix in routine ludcmp\n");
    vv[i-1]=1.0/big;
  }
  for (j=1;j<=n;j++) {
    for (i=1;i<j;i++) {
      sum=a[i-1][j-1];
      for (k=1;k<i;k++) sum -= a[i-1][k-1]*a[k-1][j-1];
      a[i-1][j-1]=sum;
    }
    big=0.0;
    for (i=j;i<=n;i++) {
      sum=a[i-1][j-1];
      for (k=1;k<j;k++)
	sum -= a[i-1][k-1]*a[k-1][j-1];
      a[i-1][j-1]=sum;
      if ( (dum=vv[i-1]*fabs(sum)) >= big) {
	big=dum;
	imax=i;
      }
    }
    if (j != imax) {
      for (k=1;k<=n;k++) {
	dum=a[imax-1][k-1];
	a[imax-1][k-1]=a[j-1][k-1];
	a[j-1][k-1]=dum;
      }
      *d = -(*d);
      vv[imax-1]=vv[j-1];
    }
    indx[j-1]=imax;
    if (a[j-1][j-1] == 0.0) a[j-1][j-1]=TINY;
    if (j != n) {
      dum=1.0/(a[j-1][j-1]);
      for (i=j+1;i<=n;i++) a[i-1][j-1] *= dum;
    }
  }
}

void lubksb(double a[equa][equa], int n, int indx[equa], double b[equa]) {

  int i,ii=0,ip,j;
  double sum;

  for (i=1;i<=n;i++) {
    ip=indx[i-1];
    sum=b[ip-1];
    b[ip-1]=b[i-1];
    if (ii)
      for (j=ii;j<=i-1;j++) sum -= a[i-1][j-1]*b[j-1];
    else if (sum) ii=i;
    b[i-1]=sum;
  }
  for (i=n;i>=1;i--) {
    sum=b[i-1];
    for (j=i+1;j<=n;j++) sum -= a[i-1][j-1]*b[j-1];
    b[i-1]=sum/a[i-1][i-1];
  }
}

/****************************************************************************
  The next two routines are a Gauss-Seidel lin. eqs solver.
*****************************************************************************/

double Dist(double x[n_grids], double y[n_grids]) {

/* Distance defined by infinity norm */

  int i;
  double temp=0.0;

  for(i=0; i<n_grids; i++)
    temp+=fabs(x[i]-y[i]);
  return temp;
} /* Dist */
  

void GaussSeidel(double A[n_grids][n_grids], double b[n_grids], double x[n_grids]) {

/* From [5] */

  int i, j, iter=0;
  double temp, x_old[n_grids];

  do {
    iter++;
    for(i=0; i<n_grids; i++)
      x_old[i]=x[i];
    for(i=1; i<=n_grids; i++) {
      temp=-A[i-1][i-1]*x[i-1];
      for(j=1; j<=n_grids; j++)
	temp+=A[i-1][j-1]*x[j-1];
      x[i-1]=(b[i-1]-temp)/A[i-1][i-1];
    } /* for i */
  } while ((Dist(x, x_old)>=EPS) && (iter<=MAX_ITER));
} /* GaussSeidel */


void initialize(char paramfile[25]) {

  int    i, l;
  FILE   *infile;
  char   filename[25];
  double noise_level, init_val[equa];

  init_diff_const();
  infile=fopen(paramfile, "r");
  fscanf(infile, "dt......= %lg\n", &dt);
  fscanf(infile, "L.......= %lg\n", &L);
  fscanf(infile, "prefix..= %s\n", prefix);
  fscanf(infile, "max_t...= %lg\n", &max_t);
  fscanf(infile, "pt1.....= %lg\n", &ptime1);
  fscanf(infile, "pt2.....= %lg\n", &ptime2);
  fscanf(infile, "period..= %lg\n", &period);
  for(l=0; l<equa; l++)
    fscanf(infile, "init-val= %lg\n", &init_val[i]);
  fscanf(infile, "noise...= %lg\n", &noise_level);
  fclose(infile);
  dx=L/n_grids;
  for(i=0; i<equa; i++) 
    R[i]=D[i]*dt/(dx*dx);
  sprintf(filename, "%s.in0", prefix);
  if ((infile=fopen(filename, "r"))!=NULL) {
    for(l=0; l<equa; l++) {
      sprintf(filename, "%s.in%d", prefix, l);
      infile=fopen(filename, "r");
      for(i=0; i<n_grids; i++) 
        fscanf(infile, "%lg", &u[i][l]);
      fclose(infile);
    }
  }
  else {
    mar_init(1, 1);
    for(l=0; l<equa; l++) {
      for(i=0; i<n_grids; i++)
	u[i][l]=init_val[l]*(1.0+noise_level*(2.0*RandUnit()-1.0));
    } /* for i */
  }
}
  
void RungeKutta(double x[n_grids][equa], double f[n_grids][equa]) {

/* From [3] */

  double k1[equa], k2[equa], k3[equa], k4[equa], temp[equa];
  int i, j;

  for(i=0; i<n_grids; i++) {
    for(j=0; j<equa; j++) {
      k1[j]=dt*reac(x[i], j);
      temp[j]=x[i][j]+k1[j]*0.5;
    }
    for(j=0; j<equa; j++) {
      k2[j]=dt*reac(temp, j);
      temp[j]=x[i][j]+k2[j]*0.5;
    }
    for(j=0; j<equa; j++) {
      k3[j]=dt*reac(temp, j);
      temp[j]=x[i][j]+k3[j];
    }
    for(j=0; j<equa; j++)
      k4[j]=x[i][j]+dt*reac(temp, j);
    for(j=0; j<equa; j++) 
      temp[j]=x[i][j]+k1[j]/6+k2[j]/3+k3[j]/3+k4[j]/6;
    for(j=0; j<equa; j++)
      f[i][j]=reac(temp, j);
  } /* for i */
} /* RungeKutta */

/*************************************************************************
  The Kronecker delta.
**************************************************************************/

int Kron(int i, int j) {

  if (i==j)
    return 1;
  else
    return 0;
}

/*************************************************************************
  The procedure implements Calahan's scheme for solving ODEs.
  See the comments for RungeKutta for details.
*************************************************************************/

void Calahan(double x[n_grids][equa], double f[n_grids][equa]) {

  const double a1=0.788675134595, b1=-1.15470053838;
  const double R1=0.75, R2=0.25;
  int i, j, l, l2;
  double B[equa][equa], Jac[equa][equa];
  double y1[equa], y2[equa], y2a[equa];
  double dummy;
  int index[equa];

  for(i=0; i<n_grids; i++) {
    calc_jac2(x[i], Jac);
    for(l=0; l<equa; l++)
      for(l2=0; l2<equa; l2++)
	B[l][l2]=(double)Kron(l, l2)-dt*a1*Jac[l][l2];
    ludcmp(B, equa, index, &dummy);
    for(l=0; l<equa; l++)
      y1[l]=dt*reac(x[i], l);
    lubksb(B, equa, index, y1);
    for(l=0; l<equa; l++)
      y2a[l]=x[i][l]+b1*y1[l];
    for(l=0; l<equa; l++)
      y2[l]=dt*reac(y2a, l);
    lubksb(B, equa, index, y2);
    for(l=0; l<equa; l++)
      y2a[l]=x[i][l]+R1*y1[l]+R2*y2[l];
    for(l=0; l<equa; l++) 
      y2a[l]=x[i][l]+R1*y1[l]+R2*y2[l];
    for(l=0; l<equa; l++)
      f[i][l]=reac(y2a, l);
  }
} /* Calahan */

/***************************************************************************
  The following routine contains the main loop.
****************************************************************************/

void iterate(void) {

/* From [1], [4] */

  double A[n_grids][n_grids], u_new[equa][n_grids];
  double f[n_grids][equa], alpha, beta;
  double a[n_grids], b[n_grids], c[n_grids], r[n_grids];
  int i, j, l, iter=0;
  FILE *outfile, *point_1[equa], *point_2[equa], *point_3[equa];
  char filename[25];
  double t;

  for(i=0; i<n_grids; i++)
    for(j=0; j<n_grids; j++)
      A[i][j]=0.0;
  for(i=0; i<equa; i++) {
    sprintf(filename, "%s-1.si%d", prefix, i);
    point_1[i]=fopen(filename, "w");
    sprintf(filename, "%s-2.si%d", prefix, i);
    point_2[i]=fopen(filename, "w");
    sprintf(filename, "%s-3.si%d", prefix, i);
    point_3[i]=fopen(filename, "w");
  }

  for(t=0.0; t<=max_t; t+=dt) {
    iter++;
/*    RungeKutta(u, f); */
    Calahan(u, f);
    for(l=0; l<equa; l++) {
      if (LINSOLVE==1) {
	for(i=1; i<(n_grids-1); i++) {
	  A[i][i-1]=-LAMBDA*R[l];
	  A[i][i]=1+2*LAMBDA*R[l];
	  A[i][i+1]=-LAMBDA*R[l];
	} /* for i */
	A[0][0]=1+2*LAMBDA*R[l];
	A[0][1]=-LAMBDA*R[l];
	A[0][n_grids-1]=-LAMBDA*R[l];
	A[n_grids-1][0]=-LAMBDA*R[l];
	A[n_grids-1][n_grids-2]=-LAMBDA*R[l];
	A[n_grids-1][n_grids-1]=1+2*LAMBDA*R[l];
	for(i=1; i<(n_grids-1); i++) 
	  b[i]=u[i-1][l]*(1-LAMBDA)*R[l]+u[i][l]*(1-2*R[l]*(1-LAMBDA))
	    +u[i+1][l]*(1-LAMBDA)*R[l]+dt*f[i][l];
	b[0]=(1-2*R[l]*(1-LAMBDA))*u[0][l]+u[1][l]*R[l]*(1-LAMBDA)
	  +u[n_grids-1][l]*R[l]*(1-LAMBDA)+dt*f[0][l];
	b[n_grids-1]=u[n_grids-2][l]*R[l]*(1-LAMBDA)
	  +u[n_grids-1][l]*(1-2*R[l]*(1-LAMBDA))+u[0][l]*R[l]*(1-LAMBDA)
	    +dt*f[n_grids-1][l];	
	GaussSeidel(A, b, u_new[l]);
      } else {
	for(i=1; i<(n_grids-1); i++) {
	  a[i]=-LAMBDA*R[l];
	  b[i]=1.0+2.0*LAMBDA*R[l];
	  c[i]=-LAMBDA*R[l];
      	  r[i]=u[i-1][l]*(1.0-LAMBDA)*R[l]+u[i][l]*(1.0-2.0*R[l]*(1.0-LAMBDA))
	    +u[i+1][l]*(1.0-LAMBDA)*R[l]+dt*f[i][l];
	} 
	alpha=-LAMBDA*R[l];
	beta=-LAMBDA*R[l];
	b[0]=1.0+2.0*LAMBDA*R[l];
	c[0]=-LAMBDA*R[l];
	a[0]=0.0;
	a[n_grids-1]=-LAMBDA*R[l];
	b[n_grids-1]=1.0+2.0*LAMBDA*R[l];
	c[n_grids-1]=0.0;
	r[0]=(1.0-2.0*R[l]*(1.0-LAMBDA))*u[0][l]+u[1][l]*R[l]*(1.0-LAMBDA)
	  +u[n_grids-1][l]*R[l]*(1.0-LAMBDA)+dt*f[0][l];
	r[n_grids-1]=u[n_grids-2][l]*R[l]*(1.0-LAMBDA)
	  +u[n_grids-1][l]*(1.0-2.0*R[l]*(1.0-LAMBDA))
	  +u[0][l]*R[l]*(1.0-LAMBDA)+dt*f[n_grids-1][l];	
	cyclic(a, b, c, alpha, beta, r, u_new[l], n_grids);
      }

      /* Output part */
#ifdef COMPACT
      if ((t>=ptime1) && (t<=ptime2)) {
        sprintf(filename, "%s.co%d", prefix, l, iter);
        outfile=fopen(filename, "a");
        for(i=0; i<n_grids; i++)
          fprintf(outfile, "%e\n", u_new[l][i]);
        fclose(outfile);
      } /* if */
      if (fmod(t, period)==0.0) {
	sprintf(filename, "%s.co%d.%d", prefix, l, iter);
	outfile=fopen(filename, "w");
	for(i=0; i<n_grids; i++)
	  fprintf(outfile, "%e\n", u_new[l][i]);
	fclose(outfile);
      } /* if */
#else
      if (fmod(t, period)==0.0) {
	sprintf(filename, "%s.co%d.%d", prefix, l, iter);
	outfile=fopen(filename, "w");
	for(i=0; i<n_grids; i++)
	  fprintf(outfile, "%e\n", u_new[l][i]);
	close(outfile);	
      }
#endif
      fprintf(point_1[l], "%e    %e\n", t, u_new[l][0]);
      fprintf(point_2[l], "%e    %e\n", t, u_new[l][n_grids % 3]);
      fprintf(point_3[l], "%e    %e\n", t, u_new[l][2*n_grids % 3]);

    } /* for l*/
    for(l=0; l<equa; l++)
      for(i=0; i<n_grids; i++)
	u[i][l]=u_new[l][i];
  } /* for t */
  for(i=0; i<equa; i++) {
    fclose(point_1[i]);
    fclose(point_2[i]);    
    fclose(point_2[i]);
  }
} /* iterate */

void main(int argc, char *argv[]) {

  if (argc!=2)
    fprintf(stderr, "Usage: %s filename\n", argv[0]);
  else {
    initialize(argv[1]);
    iterate();
  }
}






