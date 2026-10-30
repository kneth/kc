/************************************************************************* 
  Mixed - a code generator for kc. 

  CopyWrong 1993-1995 by 
  Keld Nielsen (kn@kin.kiku.dk)
  Department of Theoretical Chemistry
  University of Copenhagen
  Universitetsparken 5
  2100 Copenhagen
  Denmark
  
  See kc.tex for details

  Last updated: 8 Nov 1995 by KG
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
#include "matrix.h"
#include "quench.h"

void Mixed(void) {
  
  double  charge, temp;
  char    *name, *rename;
  Tree    tmp;
  int     i, j, l, k, no_eval, ref;
  int     num_of_spec, max_iter,num_of_reac;
  double  angle,eps, sum_left, sum_right, sumfw, sumrv,sumfw_im,sumrv_im;
  double  qd, fd;
  double  *am,*ph,*fu,*qu;
  double  *stconc, *reacfw_, *reacrv_, *flowfw, *flowrv;
  double  **jac_num;
  double  **QP, **INVQP;
  double  **reacfwds_, **reacrvds_;
  double complex z_, *values, **vectors;

  name=StringAlloc();
  rename=StringAlloc();

  num_of_spec=NoOfSpec()+NumOfDynVar()-NumOfConstraint(); 
  num_of_reac= NoOfReact();
  jac_num= MatrixAlloc(num_of_spec);
  QP= MatrixAlloc(num_of_spec);
  INVQP= MatrixAlloc(num_of_spec);
  am= VectorAlloc(num_of_spec);
  ph= VectorAlloc(num_of_spec);
  fu= VectorAlloc(num_of_spec);
  qu= VectorAlloc(num_of_spec);
  reacfw_= VectorAlloc(num_of_reac);
  reacrv_= VectorAlloc(num_of_reac);
  flowfw= VectorAlloc(num_of_reac);
  flowrv= VectorAlloc(num_of_reac);
  stconc= VectorAlloc(num_of_spec);
  values=ComplexVectorAlloc(num_of_spec);
  vectors=ComplexMatrixAlloc(num_of_spec);

  if (num_of_spec>num_of_reac) {
    reacfwds_= MatrixAlloc(num_of_spec);
    reacrvds_= MatrixAlloc(num_of_spec);
  } else {
    reacfwds_= MatrixAlloc(num_of_reac);
    reacrvds_= MatrixAlloc(num_of_reac);
  }

  InitCodeGenVar(num_of_spec, NumOfConstraint(),NoOfReact());
  GenerateRateExpr();
  GenerateJacobi();
  GenerateDiffReac();
  
  temp=GetConstant("mixed");
  if (GetError()==NoError) 
    no_eval= (int) temp;
  else
    no_eval= 0;

  temp=GetConstant("epsa");
  if (GetError()==NoError) 
    eps=temp;
  else
    eps=1.0e-18;
  temp=GetConstant("maxiter");
  if (GetError()==NoError)
    max_iter=(int)temp;
  else
    max_iter=30;
  temp=GetConstant("ref");
  if (GetError()==NoError)
    ref=(int)temp;
  else
    ref=1;

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

  for(i=0; i<num_of_reac; i++) {
    tmp=TreeCreate();
    TreeCpy(tmp, rfw[i]);
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
    reacfw_[i]=TreeEval(tmp);
    if (TreeGetError()==NoEval) {
     no_eval=1;
     fprintf(stderr, "WARNING: Was not able to evaluate element (%d) of foreward reaction vector.\n", i+1 );
      } /* if */
      TreeKill(tmp);
  } /* for i */
 
  for(i=0; i<num_of_reac; i++) {
    tmp=TreeCreate();
    TreeCpy(tmp, rrv[i]);
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
    reacrv_[i]=TreeEval(tmp);
    if (TreeGetError()==NoEval) {
     no_eval=1;
     fprintf(stderr, "WARNING: Was not able to evaluate element (%d) of reverse reaction vector.\n", i+1 );
      } /* if */
      TreeKill(tmp);
  } /* for i */

  for(i=0; i<num_of_reac; i++) {
    for(j=0; j<num_of_spec; j++) {
      tmp=TreeCreate();
      TreeCpy(tmp, rfwds[i][j]);
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
      reacfwds_[i][j]=TreeEval(tmp);
      if (TreeGetError()==NoEval) {
	no_eval=1;
	fprintf(stderr, "WARNING: Was not able to evaluate element (%d, %d) of flow reaction  matrix.\n", i+1 , j+1);
      } /* if */
      TreeKill(tmp);
    } /* for j */
  } /* for i */
 
  for(i=0; i<num_of_reac; i++) {
    for(j=0; j<num_of_spec; j++) {
      tmp=TreeCreate();
      TreeCpy(tmp, rrvds[i][j]);
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
      reacrvds_[i][j]=TreeEval(tmp);
      if (TreeGetError()==NoEval) {
	no_eval=1;
	fprintf(stderr, "WARNING: Was not able to evaluate element (%d, %d) of flow reaction  matrix.\n", i+1 , j+1);
      } /* if */
      TreeKill(tmp);
    } /* for j */
  } /* for i */

  switch (no_eval) {
  case 0:
    /* compute the eigenvalues and vectors and print them */
    Eigen(num_of_spec,1, jac_num, max_iter, eps,0, values, vectors);

    /* Print the jacobian matrix */
    printf("The Jacobian matrix is:\n");
    for(i=0; i<num_of_spec; i++) 
      for(j=0; j<num_of_spec; j++) 
	printf("J(%d, %d) = %10.14e\n", i+1, j+1, jac_num[i][j]);
    
    printf("\nThe eigenvalues are:\n");
    for(i=0; i<num_of_spec; i++)
      printf("  (%14.10e, %14.10e)\n", creal(values[i]), cimag(values[i]));

    printf("\nThe eigenvectors are:\n");
    for(i=0; i<num_of_spec; i++) {
      for(j=0; j<num_of_spec; j++)
	printf("  (%14.10e, %14.10e)\n", creal(vectors[i][j]), cimag(vectors[i][j]));
      printf("\n");
    } /* for i */

    /* control of the eigenvalues and vectors */
    for(i=0; i<num_of_spec; i++) {
      sum_left=  0.0; sum_right= 0.0;
      for(j=0; j<num_of_spec; j++) {
        sum_left +=  jac_num[i][j]*creal(vectors[i][j]);
      } /* for j */
      sum_right+= creal(values[i]*vectors[i][i]);
      printf("No. %d,REAL part: Left %e,  Right %e, Diff. %e\n",i,sum_left,sum_right,sum_left-sum_right);
      sum_left=  0.0; sum_right= 0.0;
      for(j=0; j<num_of_spec; j++) {
        sum_left +=  jac_num[i][j]*cimag(vectors[i][j]);
      } /* for j */
      sum_right+= cimag(values[i]*vectors[i][i]);
      printf("No. %d,IM   part: Left %e,  Right %e, Diff %e\n",
	     i, sum_left, sum_right, sum_left-sum_right);
    } /* for i */
    printf("\nThe eigenvalues \n");
    break;
  case 1:
    /* compute the eigenvalues and vectors and print them */
    Eigen(num_of_spec,ref, jac_num, max_iter, eps,1, values, vectors);

    i=0;
    for(l=1; l<=NoOfSpec(); l++) {
      GetSpecNo(l, name, &charge);
      if (IsSpecInConstraint(name, charge)==0) {
        RenameSpec(rename, name, charge);
        stconc[i]= GetBeginConc(name, charge);
        i += 1;
      } /* if */
    } /* for l  */
    i=0;
    do {
      if (cimag(values[i]) != 0.0) {
        for (j=0;j<num_of_spec;j++) {
          QP[j][i]= creal(vectors[i][j]);
          QP[j][i+1]= cimag(vectors[i][j]);
        }
        i+=2;
      } else {
        for (j=0;j<num_of_spec;j++) 
  	QP[j][i]= creal(vectors[i][j]);
        i+=1;
      }
    } while (i<num_of_spec);

     InversMatrix(num_of_spec,QP,INVQP);
  
     printf("\n");
     compamppha(num_of_spec, QP, am, ph);
     for(i=0;i<num_of_spec;i++)
        printf("Amp: %lf .  Phase: %lf\n", am[i], ph[i]);
     printf("\n");

     stopdata(num_of_spec,3,INVQP,stconc,qu,fu,&qd,&fd);
     for(i=0;i<num_of_spec;i++)
        printf("Que: %lf .  Phase: %lf\n", qu[i], fu[i]);
     printf("Dqu: %lf . Dphas: %lf\n\n",qd,fd);
     printf("\n");
      printf("\n Quenching \n");
    break;
  case 2:
    Eigen(num_of_spec,ref, jac_num, max_iter, eps,1, values, vectors);
    printf("Stationary flow\n");
    for(i=0;i<num_of_reac;i++)
      if (reacrv_[i] !=0.0) {
	printf("Reac. no. %d.  Foreward: %e.  Reverse: %e\n",i,reacfw_[i],reacrv_[i]);
      } else {
	printf("Reac. no. %d.  Foreward: %e.\n",i,reacfw_[i]);
      }


    printf("\n");
    i=0;
    do {
      if (cimag(values[i]) != 0.0) {
	printf("\n");
	printf("Eigenvector no. %d. Complex.\n",i+1);
	for(j=0;j<num_of_reac;j++) {
	  sumfw= 0.0; sumrv= 0.0;
	  sumfw_im= 0.0; sumrv_im= 0.0;
	  for(k=0;k<num_of_spec;k++) {
	    sumfw += reacfwds_[j][k]*creal(vectors[i][k]);
	    sumrv += reacrvds_[j][k]*creal(vectors[i][k]);
	    sumfw_im += reacfwds_[j][k]*cimag(vectors[i][k]);
	    sumrv_im += reacrvds_[j][k]*cimag(vectors[i][k]);
	  } /* for k */
	  z_= sumfw + sumfw_im*I;
	  angle= 57.296*carg(z_);
	  printf("Flow in reaction no. %d. FW: Ampl.: %e  , Phase: %e\n",j,cabs(z_),angle); 
	  if ((sumrv !=0.0)||(sumrv_im != 0.0)) {
	    z_= sumrv + sumrv_im*I;
	    angle= 57.296*carg(z_);
	    printf("Flow in reaction no. %d.RV: Ampl.: %e  , Phase: %e\n",j,cabs(z_),angle); 
	  }
	} /* for j */
	i+=2;
      } else {
	printf("\n");
	printf("Eigenvector no. %d. Real.\n",i+1);
	for(j=0;j<num_of_reac;j++) {
	  sumfw= 0.0; sumrv= 0.0;
	  for(k=0;k<num_of_spec;k++) {
	    sumfw += reacfwds_[j][k]*creal(vectors[i][k]);
	    sumrv += reacrvds_[j][k]*creal(vectors[i][k]);
	  } /* for k */
	  flowfw[j] = sumfw;
	  flowrv[j] = sumrv;
	  if (sumrv !=0.0) {
	    printf("Flow in reaction no. %d.: FW: %e  , RV: %e\n",j,flowfw[j],flowrv[j]); 
	  } else {
	    printf("Flow in reaction no. %d.: FW: %e\n",j,flowfw[j]); 
	  }
	} /* for j */
	i+=1;
      }
    } while (i<num_of_spec);
    break;
  } /* switch no_eval */
  
  VectorFree(num_of_spec,am);
  VectorFree(num_of_spec,ph);
  VectorFree(num_of_spec,fu);
  VectorFree(num_of_spec,qu);
  VectorFree(num_of_spec,stconc);
  VectorFree(num_of_reac,reacfw_);
  VectorFree(num_of_reac,reacrv_);
  MatrixFree(num_of_spec,jac_num);
  VectorFree(num_of_reac,flowfw);
  VectorFree(num_of_reac,flowrv);
  MatrixFree(num_of_spec,QP);
  MatrixFree(num_of_spec,INVQP);
  ComplexVectorFree(num_of_spec,values);
  ComplexMatrixFree(num_of_spec,vectors);
  if (num_of_spec>num_of_reac) {
    MatrixFree(num_of_spec,reacfwds_);
    MatrixFree(num_of_spec,reacrvds_);
  } else {
    MatrixFree(num_of_reac,reacfwds_);
    MatrixFree(num_of_reac,reacrvds_);
  }
} /* Mixed */
