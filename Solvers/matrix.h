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

extern double **A_fixed;
extern double *x_fixed, *s_fixed;
extern int    *p_fixed, n_fixed;

extern double  **MatrixAlloc(const int);
extern double   *VectorAlloc(const int);
extern int      *IntVectorAlloc(const int);
extern void      MatrixMul(const int, double **, double **, double **);
extern void      Transpose(const int, double **, double **);
extern void      MatrixFree(const int, double **);
extern void      VectorFree(const int, double *);
extern void      IntVectorFree(const int, int *);
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
