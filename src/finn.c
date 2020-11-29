/************************************************************************* 
  Finn - a code generator for kc. 

  CopyWrong 1993-1995 by 
  Kenneth Geisshirt (kneth@fatou.ruc.dk)
  Department of Life Sciences and Chemistry
  Roskilde University
  P.O. Box 260
  4000 Roskilde
  Denmark

  See kc.tex for details

  Last updated: 15 February 1995
*************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "config.h"
#include "symbmath.h"
#include "tableman.h"
#include "misc.h"
#include "codegen.h"
#include "eigen.h"
#include "complex.h"
#include "matrix.h"

void Finn(void) {

  double  charge, temp;
  char    *name, *rename;
  Tree    tmp;
  int     i, j, l, no_eval=0;
  int     num_of_spec, max_iter;
  double  eps;
  double  **jac_num;
  Complex *values, **vectors;

  name=StringAlloc();
  rename=StringAlloc();

  num_of_spec=NoOfSpec()+NumOfDynVar()-NumOfConstraint(); 
  jac_num=(double **)calloc(num_of_spec, sizeof(double *));
  for(i=0; i<num_of_spec; i++)
    jac_num[i]=(double *)calloc(num_of_spec, sizeof(double));
  values=ComplexVectorAlloc(num_of_spec);
  vectors=ComplexMatrixAlloc(num_of_spec);

  InitCodeGenVar(num_of_spec, NumOfConstraint());
  GenerateRateExpr();
  GenerateJacobi();
  
  temp=GetConstant("epsa");
  if (GetError()==NoError) 
    eps=temp;
  else
    temp=1.0e-18;
  temp=GetConstant("maxiter");
  if (GetError()==NoError)
    max_iter=(int)temp;
  else
    max_iter=30;

  for(i=0; i<num_of_spec; i++) {
    for(j=0; j<num_of_spec; j++) {
      tmp=TreeCreate();
      TreeCpy(tmp, jacobi[i][j]);
      for(l=1; l<=NoOfSpec(); l++) {
	GetSpecNo(l, name, &charge);
	if (IsSpecInConstraint(name, charge)==0) {
	  RenameSpec(rename, name, charge);
	  TreeSubstVar(tmp, rename, GetBeginConc(name, charge));
	} /* if */
      } /* for l  */
      for(l=1; l<=NumOfDynVar(); l++) {
	GetDynVarNo(l, name);
	TreeSubstVar(tmp, name, GetInitValue(name));
      } /* for l */
      jac_num[i][j]=TreeEval(tmp);
      if (TreeGetError()==NoEval) {
	no_eval=1;
	fprintf(stderr, "WARNING: Was not able to evaluate element (%d, %d) of Jacobian matrix.\n", i+1 , j+1);
      } /* if */
      TreeKill(tmp);
    } /* for j */
  } /* for i */
 
  if (no_eval==0) {
    /* Print the jacobian matrix */
    printf("The Jacobian matrix is:\n");
    for(i=0; i<num_of_spec; i++) 
      for(j=0; j<num_of_spec; j++) 
	printf("J(%d, %d) = %10.14e\n", i+1, j+1, jac_num[i][j]);
    
    /* compute the eigenvalues and vectors and print them */
    Eigen(num_of_spec, jac_num, max_iter, eps, values, vectors);
    printf("\nThe eigenvalues are:\n");
    for(i=0; i<num_of_spec; i++)
      printf("  (%14.10e, %14.10e)\n", values[i].re, values[i].im);
    printf("\nThe eigenvectors are:\n");
    for(i=0; i<num_of_spec; i++) {
      for(j=0; j<num_of_spec; j++)
	printf("  (%14.10e, %14.10e)\n", vectors[i][j].re, vectors[i][j].im);
      printf("\n");
    } /* for i */
  } /* if */
} /* Finn */
