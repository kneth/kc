/* SymbCont - a code generator for kc and Auto86. 
   CopyWrong by Kenneth Geisshirt, 1992, 1993

   See kc.tex for details
*/

#include <stdio.h>
#include <malloc.h>
#include <time.h>
#include "config.h"
#include "symbmath.h"
#include "tableman.h"

void AutoCode(FILE *code) {

  double      charge, temp, coeff;
  char        *name, *rename;
  time_t      timer;
  struct Tree *v[ExprSize+SymSize], *con[ConstrainSize], *r[ExprSize+SymSize];
  struct Tree *jacobi[ExprSize+SymSize][ExprSize+SymSize];
  struct Tree *v_temp, *tmp, *tmp2;
  int         i, j, k, react_no, finished, constraint, dyn, dyn2;

  name=malloc(sizeof(char));
  rename=malloc(sizeof(char));
  timer=time(&timer);

  /* process constraints */
  for(i=1; i<=NumOfConstraint(); i++) { 
    tmp=TreeCreate();
    con[i-1]=TreeCreate();
    tmp=GetConstraintNo(i, name, &charge);
    if (GetError()==NotFound) 
      fprintf(stderr, "Constraint number %d does not exist.\n", i);
    else {
      RenameSpec(rename, name, charge);
      for(j=1; j<i; j++) 
	con[j-1]=TreeSubstTree(con[j-1], rename, tmp);
      con[i-1]=TreeCpy(tmp);
    }; /* else */
    TreeKill(tmp);
  }; /* for i */
  fprintf(code, "      SUBROUTIME FUNC(ndim, u, icp, par, ijac, f, dfdu, dfdp)\n");
  fprintf(code, "C %s", ctime(&timer));
  i=NoOfSpec()+NumOfDynVar()-NumOfConstraint(); /* abuse of i */
  fprintf(code, "      IMPLICIIT REAL*8(a-h,o-z)\n");
  fprintf(code, "      DIMENSION f(ndim), dfdu(ndim, ndim), dfdp(ndim, 20)\n");
  for(i=1; i<=NoOfSpec(); i++) {
    GetSpecNo(i, name, &charge);
    RenameSpec(rename, name, charge);
    fprintf(code, "      REAL*8 %s\n", rename);
  }; /* for i */
  for(i=1; i<=NumOfDynVar(); i++) {
    GetDynVarNo(i, name);
    fprintf(code, "      REAL*8 %s\n", name);
  }; /* for i */
  for(i=1; i<=(NoOfReact()+NumOfExpr()); i++) {
    fprintf(code, "      REAL*8 r%d\n", i-1);
  }; /* for i */
  dyn=1;
  for(i=1; i<=NoOfSpec(); i++) {
    GetSpecNo(i, name, &charge);
    RenameSpec(rename, name, charge);
    constraint=IsSpecInConstraint(name, charge);
    if (constraint==0) {  
      fprintf(code, "      %s=u(%d);\n", rename, dyn);
      dyn++;
    }; /* if */
  }; /* for i*/
  for(i=1; i<=NumOfDynVar(); i++) {
    GetDynVarNo(i, name);
    fprintf(code, "      %s=par(%d)\n", name, i);
  }; /* for i */
  for(i=1; i<=NoOfSpec(); i++) {
    GetSpecNo(i, name, &charge);
    RenameSpec(rename, name, charge);
    constraint=IsSpecInConstraint(name, charge);
    if (constraint>0) {  
      fprintf(code, "      %s=", rename);
      TreePrint(con[constraint-1], 1, code);
      fprintf(code, "\n");
    }; /* if */
  }; /* for i */
  for(i=1; i<=NoOfReact(); i++) {
    r[i-1]=TreeCreate();
    react_no=GetReactNo(i-1); 
    switch (GetReactKind(react_no)) {
      case uni: if (GetRateKind(react_no, uni, 1)==1) {  
		r[i-1]=TreeAssignConst(r[i-1], GetRateConst(react_no, uni, 1));
                finished=GetFirstSpecA(react_no, name, &charge, &coeff);
		while (finished==1) {
		  if (coeff>0.0) {
		    RenameSpec(rename, name, charge);
		    temp=GetPowConstInReact(react_no, name, charge);
		    tmp=TreeCreate();
		    constraint=IsSpecInConstraint(name, charge);
		    if (constraint>0)
		      tmp=GetConstraintNo(constraint, name, &charge);
                    else
		      tmp=TreeAssignVar(tmp, rename);
		    tmp=TreePow(tmp, temp);
		    r[i-1]=TreeMul(r[i-1], tmp);
		    TreeKill(tmp);
                  }; /* if */
		  finished=GetNextSpecA(name, &charge, &coeff);
                }; /* while */
		} /* if */
		else  
		  r[i-1]=GetRateExpr(react_no, uni, 1);
		fprintf(code, "      r%d=", i-1);
		TreePrint(r[i-1], 1, code);
		fprintf(code, "\n");
                break; 
      case bi:  if (GetRateKind(react_no, bi, 1)==1) { 
		r[i-1]=TreeAssignConst(r[i-1], GetRateConst(react_no, bi, 1));
                finished=GetFirstSpecA(react_no, name, &charge, &coeff);
		while (finished==1) {
		  if (coeff>0.0) {
		    RenameSpec(rename, name, charge);
		    temp=GetPowConstInReact(react_no, name, charge);
		    tmp=TreeCreate();
		    constraint=IsSpecInConstraint(name, charge);
		    if (constraint>0)
		      tmp=GetConstraintNo(constraint, name, &charge);
                    else
		      tmp=TreeAssignVar(tmp, rename);
		    tmp=TreePow(tmp, temp);
		    r[i-1]=TreeMul(r[i-1], tmp);
		    TreeKill(tmp);
                  }; /* if */
		  finished=GetNextSpecA(name, &charge, &coeff);
                }; /* while */
		} /* if */
		else
		  r[i-1]=GetRateExpr(react_no, bi, 1);
                v_temp=TreeCreate();
                if (GetRateKind(react_no, bi, 2)==1) {
                v_temp=TreeAssignConst(v_temp, GetRateConst(react_no, bi, 2));
		finished=GetFirstSpecA(react_no, name, &charge, &coeff);
		while (finished==1) {
		  if (coeff<0.0) {
		    RenameSpec(rename, name, charge);
		    temp=GetPowConstInReact(react_no, name, charge);
		    tmp=TreeCreate();
		    constraint=IsSpecInConstraint(name, charge);
		    if (constraint>0) 
		      tmp=GetConstraintNo(constraint, name, &charge);
                    else
		      tmp=TreeAssignVar(tmp, rename);
		    tmp=TreePow(tmp, temp);
		    r[i-1]=TreeMul(r[i-1], tmp);
		    TreeKill(tmp);
                  }; /* if */
                  finished=GetNextSpecA(name, &charge, &coeff);
                }; /* while */
		} /* if */
		else 
		  v_temp=GetRateExpr(react_no, bi, 2);
		r[i-1]=TreeSub(r[i-1], v_temp);
		TreeKill(v_temp);
		fprintf(code, "      r%d=", i-1);
		TreePrint(r[i-1], 1, code);
		fprintf(code, "\n");
		break;
     case equi: fprintf(stderr, "Please use the construction: [J] = expr instead of equilibriums.\n");
		break;
   }; /* switch */
 }; /* for i */
 dyn=0;
 for(i=1; i<=NoOfSpec(); i++) {
   GetSpecNo(i, name, &charge);
   if (IsSpecInConstraint(name, charge)==0) {
     v[dyn]=TreeCreate(); 
     for(j=1; j<=NoOfReact(); j++) { 
       if (IsSpecInReact(GetReactNo(j-1), name, charge, &coeff)==1) {
	 if (j==1) {
	   fprintf(code, "      f(%d)=", dyn);
	   tmp=TreeCreate();
	   tmp=TreeAssignConst(tmp, 0);
         }; /* if */
	 v_temp=TreeCreate();
	 v_temp=TreeAssignConst(v_temp, coeff);
	 v_temp=TreeMul(v_temp, TreeCpy(r[j-1]));
	 tmp=TreeAdd(tmp, v_temp);
	 TreeKill(v_temp);
       }; /* if */
       if (j==NoOfReact()) {
         TreePrint(tmp, 1, code);
         fprintf(code, "\n");
         v[dyn]=TreeCpy(tmp);
         TreeKill(tmp);
       }; /* if */
     }; /* for j */
     dyn++;
   }; /* if */
 }; /* for i */
 for(i=1; i<=NumOfExpr(); i++) {
   tmp=TreeCreate();
   tmp=GetExprNo(i, name);
   v[dyn]=TreeCpy(tmp);
   fprintf(code, "      f(%d)=", dyn+1);
   TreePrint(tmp, 1, code);
   fprintf(code, "\n");
   dyn++;
 }; /* for i */
 fprintf(code, "      IF (ijac.EQ.0) RETURN\n");
 for(i=1; i<=NoOfReact(); i++) {
   dyn=1;
   for(j=1; j<=NoOfSpec(); j++) {
     jacobi[i-1][dyn-1]=TreeCreate(); 
     tmp=TreeCreate();
     tmp=TreeCpy(v[dyn-1]);
     GetSpecNo(j, name, &charge);
     RenameSpec(rename, name, charge);
     if (IsSpecInConstraint(name, charge)==0) {
       tmp=TreeDerive(tmp, rename);
       fprintf(code, "      dfdu(%d, %d)=", i, dyn);
       TreePrint(tmp, 1, code);
       fprintf(code, "\n");
       jacobi[i-1][dyn-1]=TreeCpy(tmp);
       dyn++;
     }; /* if */
     TreeKill(tmp);
   }; /* for j */
   dyn2=1;
   for(j=1; j<=NumOfDynVar(); j++) {
     jacobi[i-1][j-1+NoOfSpec()-NumOfConstraint()]=TreeCreate(); 
     tmp=TreeCreate();
     tmp=TreeCpy(v[dyn-1]);
     GetDynVarNo(j, name);
     if (IsVarParameter(name)==0) { 
       tmp=TreeDerive(tmp, name);
       fprintf(code, "      dfdp(%d, %d)=", i, dyn2);
       TreePrint(tmp, 1, code);
       fprintf(code, "\n");
       jacobi[i-1][dyn2-1+NoOfSpec()-NumOfConstraint()]=TreeCpy(tmp);
       dyn2++;
     }; /* if */
     TreeKill(tmp);
   }; /* for j */
 }; /* for i */
 for(i=1; i<=NumOfExpr(); i++) {
   dyn2=1;
   for(j=1; j<=NumOfExpr(); j++) {
     /* JACOBI is missing!!! */
     tmp=TreeCreate();
     tmp=TreeCpy(v[dyn-1]);
     tmp2=TreeCreate();
     tmp2=GetExprNo(j, name);
     if (IsVarParameter(name)!=0) {
       tmp=TreeDerive(tmp, name);
       fprintf(code, "      dfdx(%d, %d)=", i, dyn2);
       TreePrint(tmp, 1, code);
       fprintf(code, "\n");
       dyn2++;
     }; /* if */
     TreeKill(tmp);
     TreeKill(tmp2);
   }; /* for j */
   dyn2=1;
   for(j=1; j<=NumOfDynVar(); j++) {
     jacobi[dyn+i-1][dyn2+j-1+NoOfSpec()-NumOfConstraint()]=TreeCreate(); 
     tmp=TreeCreate();
     tmp=TreeCpy(v[dyn-1]);
     GetDynVarNo(j, name);
     if (IsVarParameter(name)==0) { 
       tmp=TreeDerive(tmp, name);
       fprintf(code, "      dfdp(%d, %d)=", i, dyn2);
       TreePrint(tmp, 1, code);
       fprintf(code, "\n");
       jacobi[dyn+i-1][dyn2-1+NoOfSpec()-NumOfConstraint()]=TreeCpy(tmp);
       dyn2++;
     }; /* if */
     TreeKill(tmp);
   }; /* for j */
 }; /* for i */
 fprintf(code, "      RETURN\n      END\n");
} /* AutoCode */ 

