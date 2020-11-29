/**************************************************************************
  A small test of matrix library.

  Kenneth Geisshirt, 29 Sep 1994
***************************************************************************/

#include <stdio.h>
#include "matrix.h"

void main(void) {

  double **A;
  double *vec, *x;
  int *index, i;

  A=MatrixAlloc(3);
  printf("A allocated\n");

  vec=VectorAlloc(3);
  printf("vec allocated\n");

  x=VectorAlloc(3);
  printf("x allocated\n");

  index=IntVectorAlloc(3);
  printf("index allocated\n");

  A[0][0]=-1.0;  A[0][1]=1.0;  A[0][2]=-4.0;
  A[1][0]=2.0;   A[1][1]=2.0;  A[1][2]=0.0;
  A[2][0]=3.0;   A[2][1]=3.0;  A[2][2]=2.0;
  
  vec[0]=0.0;  vec[1]=1.0;  vec[2]=0.5;

  LUfact(3, A, index);
  printf("LUfact\n");

  LUsubst(3, A, index, vec);
  printf("LUsubst\n");

  printf("vec = (%e, %e, %e)\n", vec[0], vec[1], vec[2]);

  A[0][0]=1.0;       A[0][1]=0.5;      A[0][2]=1.0/3.0;
  A[1][0]=1.0/3.0;   A[1][1]=1.0;      A[1][2]=0.5;
  A[2][0]=0.5;       A[2][1]=1.0/3.0;  A[2][2]=1.0;
  
  vec[0]=11.0/18.0;  vec[1]=11.0/18.0;  vec[2]=11.0/18.0;
  x[0]=0.0;          x[1]=0.0;          x[2]=0.0;

  GaussSeidel(3, A, vec, x, 1.0e-10, 100);
  printf("Gauss-Seidel\n");

  printf("x = (%e, %e, %e)\n", x[0], x[1], x[2]);

  A[0][0]=1.0;       A[0][1]=0.5;      A[0][2]=1.0/3.0;
  A[1][0]=1.0/3.0;   A[1][1]=1.0;      A[1][2]=0.5;
  A[2][0]=0.5;       A[2][1]=1.0/3.0;  A[2][2]=1.0;
  
  vec[0]=11.0/18.0;  vec[1]=11.0/18.0;  vec[2]=11.0/18.0;
  x[0]=0.0;          x[1]=0.0;          x[2]=0.0;

  Jacobi(3, A, vec, x, 1.0e-10, 100);
  printf("Jacobi\n");

  printf("x = (%e, %e, %e)\n", x[0], x[1], x[2]);

  A[0][0]=1.0;       A[0][1]=0.5;      A[0][2]=1.0/3.0;
  A[1][0]=1.0/3.0;   A[1][1]=1.0;      A[1][2]=0.5;
  A[2][0]=0.5;       A[2][1]=1.0/3.0;  A[2][2]=1.0;

  GSR(3, A);
  printf("GSR\n");

  for(i=0; i<3; i++) 
    printf("%e    %e     %e\n", A[i][0], A[i][1], A[i][2]);

  MatrixFree(3, A);
  printf("A freed\n");

  VectorFree(3, vec);
  printf("vec freed\n");

  VectorFree(3, x);
  printf("x freed\n");

  IntVectorFree(3, index);
  printf("index freed\n");
}
  
