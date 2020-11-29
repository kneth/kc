/* SymbCont - a code generator for kc and CONT. 
   This is a full operating generator.   
   CopyWrong by Kenneth Geisshirt, 1992, 1993

   See kc.tex for details
*/

#include <stdio.h>
#include <malloc.h>
#include <time.h>
#include "config.h"
#include "symbmath.h"
#include "tableman.h"
#include "codegen.h"

void SymbCont(FILE *code) {

  double charge, temp, coeff;
  char   *name, *rename;
  time_t timer;
  Tree   v[ExprSize+SymSize], con[ConstrainSize], r[ExprSize+SymSize];
  Tree   jacobi[ExprSize+SymSize+MaxParameter][ExprSize+SymSize+MaxParameter];
  Tree   v_temp, tmp, temp_tree, tmp2;
  int    i, j, k, react_no, finished, constraint, dyn, dyn2, dyn3, form;

  name=malloc(sizeof(char));
  rename=malloc(sizeof(char));
  timer=time(&timer);

  fprintf(code, "      SUBROUTINE MODEL(ndim, nvar, n, t, yy, f, g, h)\n");
  fprintf(code, "C %s", ctime(&timer));
  fprintf(code, "      IMPLICIT REAL*8(a-h,o-z)\n");
  fprintf(code, "      DIMENSION yy(ndim), f(ndim), g(ndim, nvar), h(ndim, nvar, ndim)\n");  
  i=NoOfSpec()+NumOfDynVar()-NumOfConstraint(); /* abuse of i */
  for(i=1; i<=NoOfSpec(); i++) {
    GetSpecNo(i, name, &charge);
    RenameSpec(rename, name, charge);
    fprintf(code, "    REAL*8 %s\n", rename);
  }; /* for i */
  for(i=1; i<=NumOfDynVar(); i++) {
    GetDynVarNo(i, name);
    fprintf(code, "    REAL*8 %s\n", name);
  }; /* for i */
  for(i=1; i<=(NoOfReact()+NumOfExpr()); i++) {
    fprintf(code, "    REAL*8 r%d\n", i-1);
  }; /* for i */
  fprintf(code, "      COMMOM /fixp/ par(20)\n");
  fprintf(code, "      COMMOM /varp/ alpha, beta, arg, per\n");
  dyn=1;
  for(i=1; i<=NoOfSpec(); i++) {
    GetSpecNo(i, name, &charge);
    RenameSpec(rename, name, charge);
    constraint=IsSpecInConstraint(name, charge);
    if (constraint==0) {  
      fprintf(code, "  %s=yy(%d)\n", rename, dyn);
      dyn++;
    }; /* if */
  }; /* for i*/
  for(i=1; i<=NumOfDynVar(); i++) {
    GetDynVarNo(i, name);
    fprintf(code, "  %s=par(%d)\n", name, i);
  }; /* for i */
  for(i=1; i<=NoOfSpec(); i++) {
    GetSpecNo(i, name, &charge);
    RenameSpec(rename, name, charge);
    constraint=IsSpecInConstraint(name, charge);
    if (constraint>0) {  
      fprintf(code, "  %s=", rename);
      TreePrint(con[constraint-1], 1, code);
      fprintf(code, "\n");
    }; /* if */
  }; /* for i */
  dyn=0;
  for(i=1; i<=NoOfSpec(); i++) {
    GetSpecNo(i, name, &charge);
    if (IsSpecInConstraint(name, charge)==0) {
      fprintf(code, "       f(%d)=", dyn+1); 
      TreePrint(v[dyn], 1, code);
      dyn++;
    }; /* if */
  }; /* for i */
  fprintf(code, "      IF (n.EQ.ndim) RETURN\n");
  dyn=0;
  for(i=1; i<=NoOfSpec(); i++) {
    GetSpecNo(i, name, &charge);
    if (IsSpecInConstraint(name, charge)==0) {
      dyn++;
      dyn2=1;
      for(j=1; j<=NoOfSpec(); j++) {
	tmp=TreeCreate();
	TreeCpy(tmp, v[dyn-1]);
	GetSpecNo(j, name, &charge);
	RenameSpec(rename, name, charge);
	if (IsSpecInConstraint(name, charge)==0) {
	  temp=TreeEval(jacobi[dyn-1][dyn2-1]);
	  if ((TreeGetError()==NoEval) || (temp!=0.0)) {
            fprintf(code, "       g(%d, %d)=", dyn, dyn2);
            TreePrint(tmp, 1, code);
            fprintf(code, "\n");
          } /* if */ else
            fprintf(code, "      g(%d, %d)=%e\n", dyn, dyn2, temp);
	}
      }
      for(j=1; j<=NumOfParameter(); j++) {
        tmp=TreeCreate();
        TreeCpy(tmp, v[dyn-1]);
        GetParamNo(j, name, &charge, &form);
	if (form==2)
	  RenameSpec(rename, name, charge);
	else
	  strcpy(rename, name);
	jacobi[dyn-1][NoOfSpec()-NumOfConstraint()+NumOfDynVar()+j-1]=TreeCreate();
	temp_tree=TreeCreate();
	TreeDerive(temp_tree, tmp, rename);
	TreeKill(tmp);
	tmp=TreeCreate();
	TreeCpy(tmp, temp_tree);
	TreeKill(temp_tree);
	temp=TreeEval(tmp);
	if ((TreeGetError()==NoEval) || (temp!=0.0)) {
	  fprintf(code, "       g(%d, %d)=", dyn, NoOfSpec()-NumOfConstraint()+NumOfDynVar()+j-1);
	  TreePrint(tmp, 1, code);
	  fprintf(code, "\n");
	} /* if */ else
	  fprintf(code, "      g(%d, %d)=%e\n", dyn, NoOfSpec()-NumOfConstraint()+NumOfDynVar()+j-1, temp);
	TreeCpy(jacobi[dyn-1][NoOfSpec()-NumOfConstraint()+NumOfDynVar()+j-1], tmp);
	TreeKill(tmp);
      }; /* for j */
      for(j=1; j<=NumOfDynVar(); j++) {
	jacobi[dyn-1][j-1+NoOfSpec()-NumOfConstraint()]=TreeCreate(); 
	tmp=TreeCreate();
	TreeCpy(tmp, v[dyn-1]);
	GetDynVarNo(j, name);
	temp_tree=TreeCreate();
	TreeDerive(temp_tree, tmp, rename);
	TreeKill(tmp);
	tmp=TreeCreate();
	TreeCpy(tmp, temp_tree);
	TreeKill(temp_tree);
	temp=TreeEval(tmp);
	if ((TreeGetError()==NoEval) || (temp!=0.0)) {
	  fprintf(code, "      g(%d, %d)=", dyn, j+(NoOfSpec()-NumOfConstraint()));
	  TreePrint(tmp, 1, code);
	  fprintf(code, "\n");
	}; /* if */
	TreeCpy(jacobi[dyn-1][j-1+NoOfSpec()-NumOfConstraint()], tmp);
	TreeKill(tmp);
      }; /* for j */
    }; /* if */
  }; /* for i */
  for(i=1; i<=NumOfParameter(); i++) {
    GetParamNo(i, name, &charge, &form);
    if (form==1)
      strcpy(rename, name);
    else
      RenameSpec(rename, name, charge);
    dyn=0;
    for(j=1; j<=NoOfSpec(); j++) {
      GetSpecNo(j, name, &charge);
      if (IsSpecInConstraint(name, charge)==0) {
        tmp=TreeCreate();
        TreeDerive(tmp, v[dyn], rename);
        dyn++;
        fprintf(code, "       g(%d, %d)=", dyn, NoOfSpec()-NumOfConstraint()+i);
        TreePrint(tmp, 1, code);
        fprintf(code, "\n");
        TreeKill(tmp);
      } /* if */
    } /* for j */
  } /* for i */

#ifdef HESSIAN
  fprintf(code, "      IF (n.LE.ndim*(nvar+1)) RETURN\n");
  dyn=0;
  for(i=1; i<=NoOfSpec(); i++) {
    GetSpecNo(i, name, &charge);
    if (IsSpecInConstraint(name, charge)==0) {
      dyn++;
      dyn2=1;
      for(j=1; j<=NoOfSpec(); j++) {
	tmp=TreeCreate();
	GetSpecNo(j, name, &charge);
	RenameSpec(rename, name, charge);
	if (IsSpecInConstraint(name, charge)==0) {
	  TreeCpy(tmp, jacobi[dyn-1][dyn2-1]); 
	  temp_tree=TreeCreate(); 
	  TreeDerive(temp_tree, tmp, rename);
	  dyn3=1;
	  for(k=1; k<=NoOfSpec(); k++) {
	    GetSpecNo(k, name, &charge);
	    RenameSpec(rename, name, charge);
	    if (IsSpecInConstraint(name, charge)==0) {
	      tmp2=TreeCreate();
	      TreeDerive(tmp2, temp_tree, rename);
	      fprintf(code, "      h(%d, %d, %d)=", dyn, dyn2, dyn3);
	      TreePrint(tmp2, 1, code);
	      fprintf(code, "\n");
	      dyn3++;
	      TreeKill(tmp2);
	    }; /* if */
	  }; /* for k */
	  for(k=1; k<=NumOfDynVar(); k++) {
	    GetDynVarNo(k, name);
	    tmp2=TreeCreate();
	    TreeDerive(tmp2, temp_tree, name);
	    fprintf(code, "      h(%d, %d, %d)=", dyn, dyn2, dyn3);
	    TreePrint(tmp2, 1, code);
	    fprintf(code, "\n");
	    TreeKill(tmp2);
	  }; /* for k */
	  dyn2++;
	}; /* if */
	TreeKill(tmp);
	TreeKill(temp_tree);
      }; /* for j */
      for(j=1; j<=NumOfDynVar(); j++) {
	tmp=TreeCreate();
	TreeCpy(tmp, jacobi[dyn-1][dyn2-1]);
	GetDynVarNo(j, name);
	temp_tree=TreeCreate();
	TreeDerive(temp_tree, tmp, name);
	dyn3=1;
	for(k=1; k<=NoOfSpec(); k++) {
	  GetSpecNo(k, name, &charge);
	  RenameSpec(rename, name, charge);
	  if (IsSpecInConstraint(name, charge)==0) {
	    tmp2=TreeCreate();
	    TreeDerive(tmp2, temp_tree, rename);
	    fprintf(code, "      h(%d, %d, %d)=", dyn, j+NoOfSpec()-NumOfConstraint(), dyn3);
	    TreePrint(tmp2, 1, code);
	    fprintf(code, "\n");
	    dyn3++;
	    TreeKill(tmp2);
	  }; /* if */
	}; /* for k */
	dyn2++;
	TreeKill(tmp);
	TreeKill(temp_tree);
      }; /* for j */
    }; /* if */
  }; /* for i */
#endif
  fprintf(code, "      RETURN\n      END\n");
} /* SymbCont */
