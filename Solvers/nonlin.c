/********************************************************************
  (C) Kenneth Geisshirt (kneth@fatou.ruc.dk)
      Department of Life Sciences and Chemistry
      Roskilde University
      P.O. Box 260
      4000 Roskilde
      Denmark
   
      Voice: +45 46 75 77 11, ext. 27 85
      Fax:   +45 46 75 77 21
      
  NonLin is a library containg a number of routines for solving
  nonlinear algebraic equations.
  
  
  References:
  [1]  Numerical Recipes in C, 2nd edition,
       W.H. Press, S.A. Teubolsky, W.T. Vitterling, and B.P. Flanney.
       Cambridge University Press, 1992.
  [2]  Numerical Analysis,
       D. Kincaid and W. Cheney,
       Brooks/Cole Publishing Company, 1991.
  [3]  Introduction to Numerical Analysis,
       J. Stoer and R. Bulirsch,
       Springer-Verlag, 1979.
   
  Last updated: 28 September 1995
**********************************************************************/

#include "matrix.h"
#include <stdio.h>



/*********************************************************************
  NewtonRaphson1D is the Newton-Raphson method for functions on the
  real axis, [1, pp. 362-367].
  
  Parameters:
    x0        - initial guess of the solution
    f         - the function 
    df        - the derivative of the function
    eps       - accuratecy of the solution
    max_iter  - the maximal number of iterations
**********************************************************************/

double NewtonRaphson1D(double x0, double (*f)(double), double (*df)(double),
		       double eps, int max_iter) {
  
  int     iter;      /* number of iterations */
  double  x_new, x_old;
  double  norm;
  
  x_new=x_old=x0;
  iter=0;
  do {
    iter++;
    x_old=x_new;
    x_new-=f(x_old)/df(x_old);
    norm=fabs(x_new-x_old);
  } while ((norm>=eps) || (iter<max_iter));
  return x_new;
} /* NewtonRaphson1D */


/************************************************************************
  NewtonRaphson is the Newton-Raphson method for N nonlinear equations 
  with N unknowns, [1, pp. 362-367], [3, chapter5].
*************************************************************************/

void NewtonRaphson(int N, double *x0, double *res, 
		   void (*f)(double *, double *), 
		   void (*df)(double *, double **), double eps, int max_iter) {
  
  int     iter, i;
  int     *perm;
  double  *b;
  double  **jacobi;
  
  jacobi=MatrixAlloc(N);
  perm=IntVectorAlloc(N);
  b=VectorAlloc(N);
  
  for(i=0; i<N; i++)
    res[i]=x0[i];
  
  iter=0;
  do {
    iter++;
    f(res, b);
    df(res, jacobi);
    LUfact(N, jacobi, perm);
    LUsubst(N, jacobi, perm, b);
    for(i=0; i<N; i++) 
      res[i]-=b[i];
  } while ((L2VectorNorm(N, res)>eps) && (iter<max_iter));

  IntVectorFree(N, perm);
  MatrixFree(N, jacobi);
  VectorFree(N, b);
} /* NewtonRaphson */
