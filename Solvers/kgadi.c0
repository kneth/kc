/*****************************************************************************
    This is a solver of equations like

    dc           /  d^2c    d^2c  \
    -- = f(c) + D|  ---- +  ----  |
    dt           \  dx^2    dy^2  /

    where the boundary conditions are either periodic or no-flux.

    An Altenating Direction Implicit method is used to integrate
    the spatial variables, [1] p. 855, and a 4th order Runge-Kutta,
    [1] p. 711, or Calahan's method [4] is used to integrate the 
    reaction terms (f), [1] p. 711. The linear equations are solved
    by a cyclic tridiagonal solver, [1] p. 74 or by Gauss-Seidel, [5]
    p. 356. The linear equations in Calahan are solved by a simple 
    LU-decomposition, [1] p. 44. 

    The different solvers are specified by lines in the input file macros. The
    line ODEsolve detemines which ODE solver to be used:
      1 - 4th order Runge-Kutta.
      2 - Calahan's method.
    The line LINsolve determines which solver to use in the ADI part:
      1 - Cyclic tridiagonal.
      2 - Gauss-Seidel.
    The boundary conditions are specified in the input by the line BOUNDARY, 
    which can take the following values:
      1 - Periodic. 
      2 - No flux.

    I have assumed throughout the program, that n_grids == m_grids, i.e. 
    a square lattice.

    Compilation is done by the command (on HP-UX): 
      cc -Aa +O3 -o KGadi KGadi.c model.c -lm
    Note that the +O3 option can differ from system to system. But option
    -O should be usable on all systems.

    CopyWrong, November/December 1993 by 
    Kenneth Geisshirt (kneth@osc.kiku.dk)
    Department of Theoretical Chemistry
    H. C. Orsted Institute
    Universitetsparken 5
    2100 Copenhagen
    
    Literature:

    [1] Numerical Recipes in C, 2. edition
        W. H. Press et al.
    [2] User's Manual to Kinetic Compiler v0.25
        K. Geisshirt.
    [3] Programmer's Reference Manual to Kinetic Compiler v0.25
        K. Geisshirt.
    [4] Proc. IEEE, April 1968, p. 744
        D. A. Calahan.
    [5] Matrix Computations, 1. edition
        G. H. Golub et al.
*****************************************************************************/

#include <stdio.h>
#include "model.h" /* the model is defined in this file, see [2-3]. */

#define MAX_ITER 150
#define EPS 1e-10

double u[n_grids][m_grids][equa];  /* concentrations */
double max_t, dt, ptime1, ptime2, period; /* time variables */
double R[equa];
double dx, L; /* space variable */
char   prefix[25];
double f[n_grids][n_grids][equa]; /* res. from RK */
double u_new[n_grids][n_grids][equa]; /* temp. results */
double A[n_grids][n_grids];

int ODEsolver, LINsolver, BOUNDARY; /* modes */

/*****************************************************************************
Index implements a cyclic structure.
*****************************************************************************/

int Index(int ind) {

  if (ind<0) return (ind+n_grids);
  else if (ind>=n_grids) return (ind-n_grids);
  else return ind;
} /* Index */


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

/**************************************************************************
  L_inf distance is implemented by the next function.
**************************************************************************/

double Dist(double x[n_grids], double y[n_grids]) {

  int i;
  double temp=0.0;

  for(i=0; i<n_grids; i++)
    temp+=fabs(x[i]-y[i]);
  return temp;
} /* Dist */


/**************************************************************************
  GaussSeidel implements the Gauss-Seidel scheme for solving linear 
  equations, [5].
**************************************************************************/

void GaussSeidel(double b[n_grids], double x[n_grids]) {

  int i, j, iter;
  double temp, x_old[n_grids];

  iter=0;
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



/**************************************************************************
   RungeKutta implements a 4th order Runge-Kutta method. The routines takes
   x as input. The variable x is the concentrations in each point in the
   mesh. The new (estimated) concentrations are stored in the variable
   f. The parameter h is the time step.
***************************************************************************/

void RungeKutta(double x[n_grids][m_grids][equa], 
		double f[n_grids][m_grids][equa] , double h) {

  double k1[equa], k2[equa], k3[equa], k4[equa];
  double temp[equa], temp1[equa], temp2[equa], temp3[equa];
  int i, j, l;

  for(i=0; i<n_grids; i++) {
    for(j=0; j<m_grids; j++) {
      for(l=0; l<equa; l++) {
	k1[l]=h*reac(x[i][j], l);
      temp1[l]=x[i][j][l]+k1[l]*0.5;
    }
    for(l=0; l<equa; l++) {
      k2[l]=h*reac(temp1, l);
      temp2[l]=x[i][j][l]+k2[l]*0.5;
    }
    for(l=0; l<equa; l++) {
      k3[l]=h*reac(temp2, l);
      temp3[l]=x[i][j][l]+k3[l];
    }
    for(l=0; l<equa; l++)
      k4[l]=x[i][j][l]+dt*reac(temp3, l);
    for(l=0; l<equa; l++)
      temp[l]=x[i][j][l]+k1[l]/6+k2[l]/3+k3[l]/3+k4[l]/6;
    for(l=0; l<equa; l++)
      f[i][j][l]=reac(temp, l);
    } /* for j */
  } /* for i */
} /* RungeKutta */


/************************************************************************ 
  The following function calculates the Kronecker delta.
************************************************************************/

int Kron(int i, int j) {

  if (i==j) return 1;
  else return 0;
} 

/*************************************************************************
  The procedure implements Calahan's scheme for solving ODEs. 
  See the comments for RungeKutta for details.
*************************************************************************/

void Calahan(double x[n_grids][n_grids][equa],
	     double f[n_grids][n_grids][equa], double T) {

  const double a1=0.788675134595, b1=-1.15470053838;
  const double R1=0.75, R2=0.25;
  int i, j, l, l2;
  double B[equa][equa], Jac[equa][equa];
  double y1[equa], y2[equa], y2a[equa];
  double dummy;
  int index[equa];
  
  for(i=0; i<n_grids; i++) {
    for(j=0; j<n_grids; j++) {
      calc_jac2(x[i][j], Jac);
      for(l=0; l<equa; l++)
	for(l2=0; l2<equa; l2++) 
	  B[l][l2]=(double)Kron(l, l2)-T*a1*Jac[l][l2];
      ludcmp(B, equa, index, &dummy);
      for(l=0; l<equa; l++)
	y1[l]=T*reac(x[i][j], l);
      lubksb(B, equa, index, y1);
      for(l=0; l<equa; l++) 
	y2a[l]=x[i][j][l]+b1*y1[l];
      for(l=0; l<equa; l++)
	y2[l]=T*reac(y2a, l);
      lubksb(B, equa, index, y2);
      for(l=0; l<equa; l++)
	y2a[l]=x[i][j][l]+R1*y1[l]+R2*y2[l];
      for(l=0; l<equa; l++)
	f[i][j][l]=reac(y2a, l);
    }
  }
} /* Calahan */
   
/**********************************************************************
  The procedure iterate is doing the top-level computations, i.e.
  preparing the matrices and vectors for the linear equation solvers
  and writing the data to disk.
**********************************************************************/

void iterate(void) {

  double t;             /* time */
  int    iter;          /* iteration number */
  int    i, j, l, s;       /* counters */
  char   filename[25];  /* name of a file */
  FILE   *outfile;      /* file descriptor */
  double a[n_grids], b[n_grids], c[n_grids],
         r[n_grids], alpha, beta;   /* variables for lin. eqs */
  double temp[n_grids];
  FILE   *midpoint[equa];

  iter=0;
  for(s=0; s<equa; s++) {
    sprintf(filename, "%s.si%d", prefix, s);
    midpoint[s]=fopen(filename, "w");
  }

  for(t=0.0; t<max_t; t+=dt) {
    iter++;

    if (LINsolver==2) {
      for(i=0; i<n_grids; i++)
	for(j=0; j<n_grids; j++) 
	  A[i][j]=0.0;
    }

    /* x direction */
    switch (ODEsolver) {
    case 1:
      RungeKutta(u, f, dt*0.5);  
      break;
    case 2:
      Calahan(u, f, dt*0.5); 
      break;
    }
    for(s=0; s<equa; s++) {
      for(l=0; l<n_grids; l++) {
	for(j=0; j<n_grids; j++) {
	  switch (LINsolver) {
	  case 1:
	    a[j]=-0.5*R[s];
	    b[j]=R[s]+1.0;
	    c[j]=-0.5*R[s]; 
	    break;
	  case 2:
	    A[j][Index(j-1)]=-0.5*R[s];
	    A[j][j]=1.0+R[s];
	    A[j][Index(j+1)]=-0.5*R[s];
	    break;
	  }
	}
	if (LINsolver==1) {
	  alpha=-0.5*R[s];
	  beta=-0.5*R[s];
	}
	for(j=0; j<n_grids; j++) {
          switch (BOUNDARY) {
          case 1:
	    r[j]=0.5*dt*f[j][l][s]+0.5*R[s]*u[j][Index(l+1)][s]
	      +(1.0-R[s])*u[j][Index(l)][s]+0.5*R[s]*u[j][Index(l-1)][s];
	    break;
	  case 2:
	    switch (l) {
	    case 0:
	      r[j]=0.5*dt*f[j][0][s]+0.5*R[s]*u[j][1][s]+(1.0-R[s])*u[j][0][s];
	      break;
	    case (n_grids-1):
	      r[j]=0.5*dt*f[j][n_grids-1][s]+(1.0-R[s])*u[j][n_grids-1][s]
		+0.5*R[s]*u[j][n_grids-2][s];
	      break;
	    default:
	      r[j]=0.5*dt*f[j][l][s]+0.5*R[s]*u[j][l+1][s]
		+(1.0-R[s])*u[j][l][s]+0.5*R[s]*u[j][l-1][s];
	      break;
	    } /* switch l */
	    break;
	  } /* switch BOUNDARY */
	}
	switch (LINsolver) {
	case 1:
	  if (BOUNDARY==1) 
	    cyclic(a, b, c, alpha, beta, r, temp, n_grids); 
	  else
	    tridag(a, b, c, r, temp, n_grids);
	  break;
	case 2:
	  GaussSeidel(r, temp); 
	  break;
	}
	for(j=0; j<n_grids; j++) 
	  u_new[j][l][s]=temp[j];
      } 
    }

    /* y direction */
    switch (ODEsolver) {
    case 1:
      RungeKutta(u_new, f, dt*0.5); 
      break;
    case 2:
      Calahan(u_new, f, dt*0.5); 
      break;
    }
    for(s=0; s<equa; s++) {
      for(j=0; j<n_grids; j++) {
        for(l=0; l<n_grids; l++) {
	  switch (LINsolver) {
	  case 1:
	    a[l]=-0.5*R[s];
	    b[l]=R[s]+1.0;
	    c[l]=-0.5*R[s];
	    break;
	  case 2:
	    A[l][Index(l-1)]=-0.5*R[s];
	    A[l][l]=1.0+R[s];
	    A[l][Index(l+1)]=-0.5*R[s]; 
	    break;
	  }
        }
	if (LINsolver==1) {
	  alpha=-0.5*R[s];
	  beta=-0.5*R[s];
	}
	for(l=0; l<n_grids; l++) {
	  switch (BOUNDARY) {
          case 1:
	    r[l]=0.5*dt*f[j][l][s]+0.5*R[s]*u_new[Index(j+1)][l][s]
	      +(1.0-R[s])*u_new[j][l][s]+0.5*R[s]*u_new[Index(j-1)][l][s];
	    break;
	  case 2:
	    switch (j) {
	    case 0:
	      r[l]=0.5*dt*f[0][l][s]+0.5*R[s]*u[1][l][s]+(1.0-R[s])*u[0][l][s];
	      break;
	    case (n_grids-1):
	      r[l]=0.5*dt*f[n_grids-1][l][s]+(1.0-R[s])*u[n_grids-1][l][s]
		+0.5*R[s]*u[n_grids-2][l][s];
	      break;
	    default:
	      r[l]=0.5*dt*f[j][l][s]+0.5*R[s]*u[j+1][l][s]
		+(1.0-R[s])*u[j][l][s]+0.5*R[s]*u[j-1][l][s];
	      break;
	    } /* switch l */
	    break;
	  } /* switch BOUNDARY */
	}
	switch (LINsolver) {
	case 1:
	  if (BOUNDARY==1)
	    cyclic(a, b, c, alpha, beta, r, temp, n_grids);
	  else
	    tridag(a, b, c, r, temp, n_grids);
	  break;
	case 2:
	  GaussSeidel(r, temp); 
	  break;
	}
	for(l=0; l<n_grids; l++)
	  u[j][l][s]=temp[l];
      } /* for j */
    } /* for l */

#ifdef COMPACT
    if ((t>=ptime1) && (t<=ptime2)) {
      for(l=0; l<equa; l++) {
	sprintf(filename, "%s.co%d.%d", prefix, l, iter);
	outfile=fopen(filename, "w");
	for(i=0; i<n_grids; i++)
	  for(j=0; j<m_grids; j++)
	    fprintf(outfile, "%e\n", u[i][j][l]);
	fclose(outfile);
	sprintf(filename, "gzip %s.co%d.%d", prefix, l, iter);
	system(filename);
      } /* for l */
    } /* if */
#else
    if (fmod(t, period)==0.0) {
      for(l=0; l<equa; l++) {
	sprintf(filename, "%s.co%d.%d", prefix, l, iter);
	outfile=fopen(filename, "w");
	for(i=0; i<n_grids; i++)
	  for(j=0; j<m_grids; j++)
	    fprintf(outfile, "%e\n", u[i][j][l]);
	fclose(outfile);
	sprintf(filename, "gzip %sco0%d.%d", prefix, l, iter);
	system(filename);
      } /* for l */
    }
#endif
    
    /* Saving concentrations in midpoint */
    for(s=0; s<equa; s++)
      fprintf(midpoint[s], "%e   %e\n", t, u[n_grids%2][n_grids%2][s]);
    
  } /* for t */
  
  for(s=0; s<equa; s++)
    fclose(midpoint[s]);
} /* iterate */


/********************************************************************** 
  The function initialise sets up the program. The only parameter
  is paramfile, which is the name of the setup file.
***********************************************************************/

void initialise(char paramfile[25]) {

  int i, j, l;
  FILE *infile;
  char filename[25];

  init_diff_const();
  infile=fopen(paramfile, "r");
  fscanf(infile, "dt......= %lg\n", &dt);
  fscanf(infile, "L.......= %lg\n", &L);
  fscanf(infile, "prefix..= %s\n", prefix);
  fscanf(infile, "max_t...= %lg\n", &max_t);
  fscanf(infile, "pt1.....= %lg\n", &ptime1);
  fscanf(infile, "pt2.....= %lg\n", &ptime2);
  fscanf(infile, "period..= %lg\n", &period);
  fscanf(infile, "LINsolve= %d\n", &LINsolver);
  fscanf(infile, "ODEsolve= %d\n", &ODEsolver);
  fscanf(infile, "BOUNDARY= %d\n", &BOUNDARY);
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
	for(j=0; j<m_grids; j++)
	  fscanf(infile, "%lg", &u[i][j][l]);
      fclose(infile);
    }
  }
  else {
    fprintf(stderr, "No initial conditions!\n");
    exit(1);
  }
}



/**************************************************************************
  The main program gives over the control to iterate.
**************************************************************************/

void main(int argc, char *argv[]) {

  if (argc!=2)
    fprintf(stderr, "Usage: %s filename\n", argv[0]);
  else {
    initialise(argv[1]);
    iterate();
  }
}
