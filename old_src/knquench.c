/************************************************************************* 
  KNQuench - a code generator for kc. 

  CopyWrong 1993-1995 by 
  Kenneth Geisshirt (kneth@fatou.ruc.dk)
  Department of Life Sciences and Chemistry
  Roskilde University
  P.O. Box 260
  4000 Roskilde
  Denmark

  See kc.tex for details

  Last updated: 15 May 1995 by KN
*************************************************************************/

#include <stdio.h>
#include <malloc.h>
#include <time.h>
#include "config.h"
#include "symbmath.h"
#include "tableman.h"
#include "misc.h"
#include "codegen.h"
#include "eigen.h"
#include "complex.h"
#include "matrix.h"

void KnQuench(void) {

  double  charge, temp;
  char    *name, *rename;
  Tree    tmp;
  int     i, j, l, k, no_eval=1, ref;
  int     num_of_spec, max_iter,num_of_reac;
  double  angle,eps, sum_left, sum_right, sumfw, sumrv,sumfw_im,sumrv_im;
  double  qd, fd;
  double  *am,*ph,*fu,*qu;
  double  **jac_num;
  double  **QP, **INVQP;
  double  *stconc, *reacfw_, *reacrv_, *flowfw, *flowrv;
  double  **reacfwds_, **reacrvds_;
  double  leftsc,rightsc;
  Complex z_,*values, **vectors;

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
 
  if (no_eval==1) {
    /* compute the eigenvalues and vectors and print them */
    Eigen(num_of_spec,ref, jac_num, max_iter, eps,1, values, vectors);

    /* Print the jacobian matrix */
    printf("The Jacobian matrix is:\n");
    for(i=0; i<num_of_spec; i++) 
      for(j=0; j<num_of_spec; j++) 
	printf("J(%d, %d) = %10.14e\n", i+1, j+1, jac_num[i][j]);
    
    printf("\nThe eigenvalues are:\n");
    for(i=0; i<num_of_spec; i++) {
      if (values[i].im !=0.0) {
	 printf("  (%14.10e, %14.10e)\n", values[i].re, values[i].im);
      } else {
	 printf("  (%14.10e)\n", values[i].re);
      }
      }

    printf("\nThe eigenvectors are:\n");
    for(i=0; i<num_of_spec; i++) {
      for(j=0; j<num_of_spec; j++) {
        if (values[i].im !=0.0) {
          printf("  (%14.10e, %14.10e)\n", vectors[i][j].re, vectors[i][j].im);
        } else {
          printf("  (%14.10e)\n", vectors[i][j].re);
	}
      }
      printf("\n");
    } /* for i */

  /* control of the eigenvalues and vectors */
  for(i=0; i<num_of_spec; i++) {
    sum_left=  0.0; sum_right= 0.0;
    for(j=0; j<num_of_spec; j++) {
      sum_left +=  jac_num[i][j]*vectors[i][j].re;
    } /* for j */
    sum_right+= values[i].re*vectors[i][i].re-values[i].im*vectors[i][i].im;
    printf("No. %d,REAL part: Left %e,  Right %e , Diff. %e\n",i,sum_left,sum_right,sum_left-sum_right);
    sum_left=  0.0; sum_right= 0.0;
    for(j=0; j<num_of_spec; j++) {
      sum_left +=  jac_num[i][j]*vectors[i][j].im;
    } /* for j */
    sum_right+= values[i].re*vectors[i][i].im+values[i].im*vectors[i][i].re;
    printf("No. %d,IM   part: Left %e,  Right %e , Diff %e\n",i,sum_left,sum_right,sum_left-sum_right);
    } /* for i */

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
    if (values[i].im != 0.0) {
      for (j=0;j<num_of_spec;j++) {
	QP[j][i]= vectors[i][j].re;
	QP[j][i+1]= vectors[i][j].im;
      }
      i+=2;
    } else {
      for (j=0;j<num_of_spec;j++) 
	QP[j][i]= vectors[i][j].re;
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
     if (values[i].im != 0.0) {
       printf("\n");
       printf("Eigenvector no. %d. Complex.\n",i+1);
       for(j=0;j<num_of_reac;j++) {
	 sumfw= 0.0; sumrv= 0.0;
	 sumfw_im= 0.0; sumrv_im= 0.0;
	 for(k=0;k<num_of_spec;k++) {
	   sumfw += reacfwds_[j][k]*vectors[i][k].re;
	   sumrv += reacrvds_[j][k]*vectors[i][k].re;
	   sumfw_im += reacfwds_[j][k]*vectors[i][k].im;
	   sumrv_im += reacrvds_[j][k]*vectors[i][k].im;
	 } /* for k */
	 ComplexAssign(sumfw,1.0*sumfw_im,&z_);
	 angle= 57.296*ComplexArg(z_);
	 printf("Flow in reaction no. %d. FW: Ampl.: %e  , Phase: %e\n",j,ComplexNorm(z_),angle); 
         if ((sumrv !=0.0)||(sumrv_im != 0.0)) {
	   ComplexAssign(sumrv,1.0*sumrv_im,&z_);
	 angle= 57.296*ComplexArg(z_);
	   printf("Flow in reaction no. %d.RV: Ampl.: %e  , Phase: %e\n",j,ComplexNorm(z_),angle); 
         }
       } /* for j */
       i+=2;
     } else {
       printf("\n");
       printf("Eigenvector no. %d. Real.\n",i+1);
       for(j=0;j<num_of_reac;j++) {
	 sumfw= 0.0; sumrv= 0.0;
	 for(k=0;k<num_of_spec;k++) {
	   sumfw += reacfwds_[j][k]*vectors[i][k].re;
	   sumrv += reacrvds_[j][k]*vectors[i][k].re;
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
  } /* if */

  /* stoikiometrisk matrix */
   for(i=0;i<num_of_reac;i++) {
    for(j=0; j< num_of_spec; j++) {
    printf("TEST 1\n");
      GetSpecNo(j, name, &charge);
    printf("TEST 2\n");
    printf("Reac no. %d   %s  %e\n",i,name,charge);
    printf("TEST 3\n");
      leftsc=  GetPowConstInReact(i,name,charge,0);
      /*
      leftsc=  GetCoeffInReact(i,name,charge,0);
      */
      rightsc= GetCoeffInReact(i,name,charge,1);
    printf("TEST 4\n");
      printf("Reac.no.: %d  . Left: %e  Right: %e\n",i,leftsc,rightsc);
    } 
    printf("\n");
  }

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
} /* KnQuench */
