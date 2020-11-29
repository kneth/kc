/****************************************************************************
  Misc. routines for manipulating matrices and solving linear equations.
  Matrices are assumed to be declared as **double and allocated by the
  function MatrixAlloc. A matrix can be freed by MatrixFree. Similar for
  vectors.

  CopyWrong 1994 by
  Kenneth Geisshirt (kneth@fatou.ruc.dk)
  Department of Life Sciences and Chemistry
  Roskilde University
  P.O. Box 260
  4000 Roskilde
  Denmark

  Last updated: 9 May 1995 by KN
*****************************************************************************/

#ifndef _MATRIX_LIB_
#define _MATRIX_LIB_

#include "complex.h"

#ifndef MALLOCTYPE         /* init. the type used by free(3) */
#ifdef _PLATFORM_DOS_
#define MALLOCTYPE void
#endif
#ifdef _PLATFORM_GCC_
#define MALLOCTYPE void
#endif
#ifdef _PLATFORM_HPUX_
#define MALLOCTYPE char
#endif
#ifdef _PLATFORM_CONVEX_
#define MALLOCTYPE void
#endif
#ifdef _PLATFORM_ULTRIX_
#define MALLOCTYPE char
#endif
#ifdef _PLATFORM_SGI_
#define MALLOCTYPE void
#endif
#ifdef _PLATFORM_LINUX_
#define MALLOCTYPE void
#endif
#ifdef _PLATFORM_AIX_
#define MALLOCTYPE void
#endif
#endif

double **A_fixed;
double *x_fixed, *s_fixed;
int    *p_fixed, n_fixed;

extern double  **MatrixAlloc(const int);
extern double   *VectorAlloc(const int);
extern int      *IntVectorAlloc(const int);
extern Complex  *ComplexVectorAlloc(const int);
extern Complex **ComplexMatrixAlloc(const int);
extern void      MatrixMul(const int, double **, double **, double **);
extern void      Transpose(const int, double **, double **);
extern void      MatrixFree(const int, double **);
extern void      VectorFree(const int, double *);
extern void      IntVectorFree(const int, int *);
extern void      ComplexMatrixFree(const int, Complex **);
extern void      ComplexVectorFree(const int, Complex *);
extern void      LUfact(const int, double **, int *);
extern void      LUsubst(const int, double **, int *, double *);
extern void      LUfact_fixed(void);
extern void      LUsubst_fixed(double *);
extern void      LU_fixed_init(int);
extern void      Tridiag(const int, double *, double *, double *, double *);
extern void      GaussSeidel(const int, double **, double *, double *, 
			     double, int); 
extern void      Jacobi(const int, double **, double *, double *, double, 
			int); 
extern double    DotProd(const int, double *, double *);
extern void      MatrixVecProd(const int, double **, double *, double *);
extern void      MatrixCopy(const int, double **, double **);
extern void      GSR(const int, double **);
extern double    L2VectorNorm(const int, double *);
extern void      InversMatrix(const int, double **, double **);
#endif
