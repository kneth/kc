/**************************************************************************
  A small test of the eigenvalue library. The test example comes from
  many sources.

  CopyWrong 1994-1995 by Kenneth Geisshirt.

  Last updated: 9 November 1995
***************************************************************************/

#include <stdio.h>

#include "eigen.h"
#include "matrix.h"
#include "complex.h"

void main(void) {

  double **A;                /* matrices */
  Complex *val, **vecs;
  int     i, j;

  A=MatrixAlloc(3);
  val=ComplexVectorAlloc(3);
  vecs=ComplexMatrixAlloc(3);

  printf("Matrix 1:\n");
  A[0][0]=1.000;   A[0][1]=0.000;   A[0][2]=0.010; 
  A[1][0]=0.100;   A[1][1]=1.000;   A[1][2]=0.000;
  A[2][0]=0.000;   A[2][1]=1.000;   A[2][2]=1.000;

  Eigen(3, 1, A, 150, 1.0e-15, 0, val, vecs);

  printf("Eigenvalues:\n");
  for(i=0; i<3; i++)
    printf("  (%e, %e)\n", val[i].re, val[i].im);
  printf("\nEigenvectors:\n");
  for(i=0; i<3; i++) {
    for(j=0; j<3; j++)
      printf("  (%e, %e)\n", vecs[i][j].re, vecs[i][j].im);
    printf("\n");
  }
    
  printf("Matrix 2:\n");
  A[0][0]=1.000;   A[0][1]=0.000;   A[0][2]=0.000; 
  A[1][0]=0.000;   A[1][1]=2.000;   A[1][2]=0.000;
  A[2][0]=0.000;   A[2][1]=0.000;   A[2][2]=3.000;

  Eigen(3, 1, A, 150, 1.0e-15, 0, val, vecs);

  printf("Eigenvalues:\n");
  for(i=0; i<3; i++)
    printf("  (%e, %e)\n", val[i].re, val[i].im);
  printf("\nEigenvectors:\n");
  for(i=0; i<3; i++) {
    for(j=0; j<3; j++)
      printf("  (%e, %e)\n", vecs[i][j].re, vecs[i][j].im);
    printf("\n");
  }
}

