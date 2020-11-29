/****************************************************************************
  CodeGen - code generator routines for kc.
  
  CopyWrong 1992-199 by 
  Kenneth Geisshirt (kneth@fatou.ruc.dk)        Keld Nielsen (kn@kin.kiku.dk)
  Dept. of Life Sciences and Chemistry         Dept. of Theoretical Chemistry
  Roskilde University                                University of Copenhagen
  P.O. Box 260                                           Universitetsparken 5
  DK-4000 Roskilde                                    DK-2100 Copenhagen East

  See kc.tex for details.

  Last updated: 24 May 1996 by KG
*****************************************************************************/

#include <stdio.h>
#include <string.h>
#include "tableman.h"
#include "symbmath.h"
#include "config.h"
#include "codegen.h"
#include "misc.h"

#ifdef _USE_GARBAGE_COL_
#  include <gc.h>
#else
#  include <stdlib.h>
#endif

/***************************************************************************** 
  Code generator independent routines (exported).
  All routines are using global variables defined in the header file.
*****************************************************************************/

void InitCodeGenVar(int n, int m, int r) {

  int i, j, k;

#ifdef _USE_GARBAGE_COL_
  v=(Tree *)GC_malloc(n*sizeof(Tree));
  con=(Tree *)GC_malloc(m*sizeof(Tree));
  jacobi=(Tree **)GC_malloc(n*sizeof(Tree *));
  for(i=0; i<n; i++) 
    jacobi[i]=(Tree *)GC_malloc(n*sizeof(Tree));
  hess=(Tree ***)GC_malloc(n*sizeof(Tree **));
  for(i=0; i<n; i++) {
    hess[i]=(Tree **)GC_malloc(n*sizeof(Tree *));
    for(j=0; j<n; j++)
      hess[i][j]=(Tree *)GC_malloc(n*sizeof(Tree));
  } /* for i */
  keld=(Tree ****)GC_malloc(n*sizeof(Tree ***));
  for(i=0; i<n; i++) {
    keld[i]=(Tree ***)GC_malloc(n*sizeof(Tree **));
    for(j=0; j<n; j++) {
      keld[i][j]=(Tree **)GC_malloc(n*sizeof(Tree *));
      for(k=0; k<n; k++)
	keld[i][j][k]=(Tree *)GC_malloc(n*sizeof(Tree));
    } /* for j */
  } /* for i */
#else
  v=(Tree *)calloc(n, sizeof(TreeNode));
  con=(Tree *)calloc(m, sizeof(TreeNode));
  jacobi=(Tree **)calloc(n, sizeof(TreeNode *));
  for(i=0; i<n; i++) 
    jacobi[i]=(Tree *)calloc(n, sizeof(TreeNode));
  hess=(Tree ***)calloc(n, sizeof(TreeNode **));
  for(i=0; i<n; i++) {
    hess[i]=(Tree **)calloc(n, sizeof(TreeNode *));
    for(j=0; j<n; j++)
      hess[i][j]=(Tree *)calloc(n, sizeof(TreeNode));
  } /* for i */
  keld=(Tree ****)calloc(n, sizeof(TreeNode ***));
  for(i=0; i<n; i++) {
    keld[i]=(Tree ***)calloc(n, sizeof(TreeNode **));
    for(j=0; j<n; j++) {
      keld[i][j]=(Tree **)calloc(n, sizeof(TreeNode *));
      for(k=0; k<n; k++)
	keld[i][j][k]=(Tree *)calloc(n, sizeof(TreeNode));
    } /* for j */
  } /* for i */
  rfw=(Tree *)calloc(r, sizeof(TreeNode));
  rrv=(Tree *)calloc(r, sizeof(TreeNode));
  if(r>n) {
    rfwds=(Tree **)calloc(r, sizeof(TreeNode *));
    for(i=0; i<r; i++) 
      rfwds[i]=(Tree *)calloc(r, sizeof(TreeNode));
    rrvds=(Tree **)calloc(r, sizeof(TreeNode *));
    for(i=0; i<r; i++) 
      rrvds[i]=(Tree *)calloc(r, sizeof(TreeNode));
  } else {
    rfwds=(Tree **)calloc(n, sizeof(TreeNode *));
    for(i=0; i<n; i++) 
      rfwds[i]=(Tree *)calloc(n, sizeof(TreeNode));
    rrvds=(Tree **)calloc(n, sizeof(TreeNode *));
    for(i=0; i<n; i++) 
      rrvds[i]=(Tree *)calloc(n, sizeof(TreeNode));
  }
#endif
} /* InitCodeGenVar */


void GenerateRateExpr(void) {

  double charge, temp, coeff;
  char   *name, *rename, *strtmp;
  Tree   v_temp, tmp, temp_tree, tree_temp;
  Tree   *v_store, v_tree, *r;
  int    i, j, react_no, finished, constraint, dyn, dyn2;

#ifdef _USE_GARBAGE_COL_
  v_store=(Tree *)GC_malloc((NoOfSpec()-NumOfConstraint()+NumOfDynVar()+IsNonAutoSystem())*sizeof(Tree));
  r=(Tree *)GC_malloc((ExprSize+SymSize)*sizeof(Tree));
#else
  v_store=(Tree *)calloc(NoOfSpec()-NumOfConstraint()+NumOfDynVar()+IsNonAutoSystem(), sizeof(Tree));
  r=(Tree *)calloc((ExprSize+SymSize), sizeof(Tree));
#endif
  name=StringAlloc();
  rename=StringAlloc();
  strtmp=StringAlloc();

  /* process constraints */
  for(i=1; i<=NumOfConstraint(); i++) { 
    tmp=TreeCreate();
    con[i-1]=TreeCreate();
    GetConstraintNo(i, name, &charge, tmp);
    RenameSpec(rename, name, charge);
    for(j=1; j<i; j++) 
      TreeSubstTree(con[j-1], rename, tmp);
    TreeCpy(con[i-1], tmp);
    TreeKill(tmp);
  } /* for i */
  for(i=1; i<=NoOfReact(); i++) {
    r[i-1]=TreeCreate();
    rfw[i-1]=TreeCreate();
    rrv[i-1]=TreeCreate();
    react_no=GetReactNo(i-1); 
    switch (GetReactKind(react_no)) {
    case uni: 
      if (GetRateKind(react_no, uni, 1)==1) {
	GetRateConst(react_no, uni, 1, r[i-1]);
	finished=GetFirstSpecA(react_no, name, &charge, &coeff, 0);
	while (finished==1) {
	  if (coeff!=0.0) {
	    RenameSpec(rename, name, charge);
	    temp=GetPowConstInReact(react_no, name, charge, 0);
	    temp_tree=TreeCreate();
	    TreeAssignConst(temp_tree, temp);
	    tmp=TreeCreate();
	    constraint=IsSpecInConstraint(name, charge);
	    if (constraint>0)
	      GetConstraintNo(constraint, name, &charge, tmp);
	    else
	      TreeAssignVar(tmp, rename);
	    TreePow(tmp, temp_tree);
	    TreeKill(temp_tree);
	    TreeMul(r[i-1], tmp);
	    TreeKill(tmp);
	  } /* if */
	  finished=GetNextSpecA(name, &charge, &coeff, 0);
	} /* while */
      } /* if */
      else {  
	GetRateExpr(react_no, uni, 1, r[i-1]);
	for(j=1; j<=NoOfSpec(); j++) {
	  GetSpecNo(j, name, &charge);
	  constraint=IsSpecInConstraint(name, charge);
	  if (constraint>0) {
	    RenameSpec(rename, name, charge);
	    TreeSubstTree(r[i-1], rename, con[constraint-1]);
	  } /* if */
	} /* for j */
      } /* else */
      TreeCpy(rfw[i-1],r[i-1]);
      TreeAssignConst(rrv[i-1],0.0);
      break; /* case uni */
    case bi:  
      if (GetRateKind(react_no, bi, 1)==1) { 
	GetRateConst(react_no, bi, 1, r[i-1]);
	finished=GetFirstSpecA(react_no, name, &charge, &coeff, 0);
	while (finished==1) {
	  if (coeff!=0.0) {
	    RenameSpec(rename, name, charge);
	    temp=GetPowConstInReact(react_no, name, charge, 0);
	    temp_tree=TreeCreate();
	    TreeAssignConst(temp_tree, temp);
	    tmp=TreeCreate();
	    constraint=IsSpecInConstraint(name, charge);
	    if (constraint>0)
	      GetConstraintNo(constraint, name, &charge, tmp);
	    else
	      TreeAssignVar(tmp, rename);
	    TreePow(tmp, temp_tree);
	    TreeMul(r[i-1], tmp);
	    TreeKill(tmp);
	    TreeKill(temp_tree);
	  } /* if */
	  finished=GetNextSpecA(name, &charge, &coeff, 0);
	} /* while */
      } /* if */
      else {
	GetRateExpr(react_no, bi, 1, r[i-1]);
        for(j=1; j<=NoOfSpec(); j++) {
          GetSpecNo(j, name, &charge);
          constraint=IsSpecInConstraint(name, charge);
          if (constraint>0) {
            RenameSpec(rename, name, charge);
            TreeSubstTree(r[i-1], rename, con[constraint-1]);
          } /* if */
        } /* for j */
      } /* else */
      TreeCpy(rfw[i-1],r[i-1]);
      v_temp=TreeCreate();
      if (GetRateKind(react_no, bi, 2)==1) {
	GetRateConst(react_no, bi, 2, v_temp);
	finished=GetFirstSpecA(react_no, name, &charge, &coeff, 1);
	while (finished==1) {
	  if (coeff!=0.0) {
	    RenameSpec(rename, name, charge);
	    temp=GetPowConstInReact(react_no, name, charge, 1);
	    tmp=TreeCreate();
	    temp_tree=TreeCreate();
	    TreeAssignConst(temp_tree, temp);
	    constraint=IsSpecInConstraint(name, charge);
	    if (constraint>0) 
	      GetConstraintNo(constraint, name, &charge, tmp);
	    else
	      TreeAssignVar(tmp, rename);
	    TreePow(tmp, temp_tree);
	    TreeMul(v_temp, tmp);
	    TreeKill(tmp);
	    TreeKill(temp_tree);
	  } /* if */
	  finished=GetNextSpecA(name, &charge, &coeff, 1);
	} /* while */
      } /* if */
      else {
	GetRateExpr(react_no, bi, 2, v_temp);
        for(j=1; j<=NoOfSpec(); j++) {
          GetSpecNo(j, name, &charge);
          constraint=IsSpecInConstraint(name, charge);
          if (constraint>0) {
            RenameSpec(rename, name, charge);
            TreeSubstTree(v_temp, rename, con[constraint-1]);
          } /* if */
        } /* for j */
      } /* else */
      TreeSub(r[i-1], v_temp);
      TreeCpy(rrv[i-1],v_temp);
      TreeKill(v_temp);
      break; /* case bi */
    case equi: 
      fprintf(stderr, "Please use the construction: [J] = expr instead of equilibriums.\n");
      break;
    }; /* switch */
  }; /* for i */
  dyn=0;
  for(i=1; i<=NoOfSpec(); i++) {
    GetSpecNo(i, name, &charge);
    if (IsSpecInConstraint(name, charge)==0) {
      v_store[dyn]=TreeCreate(); 
      for(j=1; j<=NoOfReact(); j++) {
	if (j==1) {
	  tmp=TreeCreate();
	  TreeAssignConst(tmp, 0.0);
	} /* if */
	if (IsSpecInReact(GetReactNo(j-1), name, charge, &coeff)==1) {
	  v_temp=TreeCreate();
	  TreeAssignConst(v_temp, -coeff);
	  TreeMul(v_temp, r[j-1]);
	  TreeAdd(tmp, v_temp);
	  TreeKill(v_temp);
	} /* if */
	if (j==NoOfReact()) {
	  TreeCpy(v_store[dyn], tmp);
	  TreeKill(tmp);
	} /* if */
      } /* for j */
      dyn++;
    } /* if */
  } /* for i */
  for(i=1; i<=NumOfDynVar(); i++) {
    v_store[i+NoOfSpec()-NumOfConstraint()-1]=TreeCreate();
    GetExprNo(i, name, v_store[i+NoOfSpec()-NumOfConstraint()-1]);
  } /* for i */
  for(i=0; i<(NoOfSpec()-NumOfConstraint()+NumOfDynVar()+IsNonAutoSystem()); 
      i++) {
    v[i]=TreeCreate();
    TreeCpy(v[i], v_store[i]);
    TreeKill(v_store[i]);
  } /* for i */

  StringFree(name);
  StringFree(rename);
} /* GenerateRateExpr */


void GenerateJacobi(void) {

  double charge, temp, coeff;
  char   *name, *rename, *strtmp;
  Tree   v_temp, tmp, temp_tree;
  int    i, j, react_no, finished, constraint, dyn, dyn2, dyn3;
  int    jj;

  name=StringAlloc();
  rename=StringAlloc();
  strtmp=StringAlloc();

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
	  jacobi[dyn-1][dyn2-1]=TreeCreate(); 
	  temp_tree=TreeCreate();
	  TreeDerive(temp_tree, tmp, rename);
	  TreeKill(tmp);
	  tmp=TreeCreate();
	  TreeCpy(tmp, temp_tree);
	  TreeKill(temp_tree);
	  temp=TreeEval(tmp);
	  TreeCpy(jacobi[dyn-1][dyn2-1], tmp);
	  dyn2++;
	} /* if */
	TreeKill(tmp);
      } /* for j */
      for(j=1; j<=NumOfDynVar(); j++) {
	jacobi[dyn-1][j-1+NoOfSpec()-NumOfConstraint()]=TreeCreate(); 
	tmp=TreeCreate();
	TreeCpy(tmp, v[dyn-1]);
	GetDynVarNo(j, name);
	temp_tree=TreeCreate();
	TreeDerive(temp_tree, tmp, name);
	TreeKill(tmp);
	tmp=TreeCreate();
	TreeCpy(tmp, temp_tree);
	TreeKill(temp_tree);
	TreeCpy(jacobi[dyn-1][j-1+NoOfSpec()-NumOfConstraint()], tmp);
	TreeKill(tmp);
      } /* for j */
    } /* if */
  } /* for i */
  for(i=1; i<=NumOfDynVar(); i++) {
    dyn=0;
    for(j=1; j<=NoOfSpec(); j++) {
      GetSpecNo(j, name, &charge);
      RenameSpec(rename, name, charge);
      if (IsSpecInConstraint(name, charge)==0) {
	dyn++;
	jacobi[i+NoOfSpec()-NumOfConstraint()-1][dyn-1]=TreeCreate();
	TreeDerive(jacobi[i+NoOfSpec()-NumOfConstraint()-1][dyn-1], 
		   v[i+NoOfSpec()-NumOfConstraint()-1], rename);
      } /* if */
    } /* for j */
    for(j=1; j<=NumOfDynVar(); j++) {
      GetDynVarNo(j, name);
      jacobi[i+NoOfSpec()-NumOfConstraint()-1][j+NoOfSpec()-NumOfConstraint()-1]=TreeCreate();
      TreeDerive(jacobi[i+NoOfSpec()-NumOfConstraint()-1][j+NoOfSpec()-NumOfConstraint()-1], v[i+NoOfSpec()-NumOfConstraint()-1], name);
    } /* for j */
  } /* for i */

  StringFree(name);
  StringFree(rename);
  StringFree(strtmp);
} /* GenerateJacobi */

void GenerateHessian(void) {

  /* It is assumed that the jacobian is generated. */

  double charge;
  char   *name, *rename;
  int    i, j, l, dyn;
  int NumbOfDynVars;

  name=StringAlloc();
  rename=StringAlloc();

  NumbOfDynVars=NoOfSpec()+NumOfDynVar()-NumOfConstraint();
  for(i=0; i<NumbOfDynVars; i++)
    for(j=0; j<NumbOfDynVars; j++) {
      dyn=0;
      for(l=1; l<=NoOfSpec(); l++) {
	GetSpecNo(l, name, &charge);
	if (IsSpecInConstraint(name, charge)==0) {
	  RenameSpec(rename, name, charge);
	  hess[i][j][dyn]=TreeCreate();
	  TreeDerive(hess[i][j][dyn], jacobi[i][j], rename);
	  dyn++;
	} /* if */
      } /* for l */
      for(l=1; l<=NumOfDynVar(); l++) {
	GetDynVarNo(l, name);
	hess[i][j][i+NoOfSpec()-NumOfConstraint()-1]=TreeCreate();
	TreeDerive(hess[i][j][i+NoOfSpec()-NumOfConstraint()-1], jacobi[i][j],
		   name);
      } /* for l */
    } /* for j */
  StringFree(name);
  StringFree(rename);
} /* GenerateHessian */

void EvalJacobian(double *x, double **jac) {

  int     dyn, dyn1, dyn2, i, j, l;
  Tree    tmp;
  double  charge;
  char    *name, *rename;

  name=StringAlloc();
  rename=StringAlloc();
  dyn=1;
  for(i=1; i<=NoOfSpec(); i++) {
    GetSpecNo(i, name, &charge);
    if (IsSpecInConstraint(name, charge)==0) {
      dyn2=1;
      for(j=1; j<=NoOfSpec(); j++) {
        GetSpecNo(j, name, &charge);
        if (IsSpecInConstraint(name, charge)==0) {
          tmp=TreeCreate();
          TreeCpy(tmp, jacobi[dyn-1][dyn2-1]);
          for(l=1; l<=NoOfSpec(); l++) {
            GetSpecNo(l, name, &charge);
            RenameSpec(rename, name, charge);
            TreeSubstVar(tmp, rename, x[l-1]);
          } /* for l  */
          for(l=1; l<=NumOfDynVar(); l++) {
            GetDynVarNo(l, name);
            TreeSubstVar(tmp, name, x[l+NoOfSpec()-NumOfConstraint()-1]);
          } /* for l */
          jac[dyn-1][dyn2-1]=TreeEval(tmp);
          if (TreeGetError()==NoEval)
            fprintf(stderr, "WARNING: Was not able to evaluate element (%d, %d) of Jacobian.\n", dyn , dyn2);
          TreeKill(tmp);
          dyn2++;
        } /* if */
      } /* for j */
      dyn++;
    } /* if */
  } /* for i */
  StringFree(name);
  StringFree(name);
} /* EvalJacobian */


void GenerateKeldian(void) {

  /* It is assumed that the hessian is generated. */

  double charge;
  char   *name, *rename;
  int    i, j, k, l, dyn;
  int    NumbOfDynVars;

  name=StringAlloc();
  rename=StringAlloc();

  NumbOfDynVars=NoOfSpec()+NumOfDynVar()-NumOfConstraint();
  for(i=0; i<NumbOfDynVars; i++)
    for(j=0; j<NumbOfDynVars; j++)
      for(k=0; k<NumbOfDynVars; k++) {
	dyn=0;
	for(l=1; l<=NoOfSpec(); l++) {
	  GetSpecNo(l, name, &charge);
	  if (IsSpecInConstraint(name, charge)==0) {
	    RenameSpec(rename, name, charge);
	    keld[i][j][k][dyn]=TreeCreate();
	    TreeDerive(keld[i][j][k][dyn], hess[i][j][k], rename);
	    dyn++;
	  } /* if */
	} /* for l */
	for(l=1; l<=NumOfDynVar(); l++) {
	  GetDynVarNo(l, name);
	  keld[i][j][k][i+NoOfSpec()-NumOfConstraint()-1]=TreeCreate();
	  TreeDerive(keld[i][j][k][i+NoOfSpec()-NumOfConstraint()-1], 
		     hess[i][j][k], name);
	} /* for l */
      } /* for k */
  StringFree(name);
  StringFree(rename);
} /* GenerateKeldian */

void GenerateDiffReac(void) {

  /* It is assumed that the vectors rfw and rrv are generated. */

  double charge;
  char   *name, *rename;
  int    i, j, NumInd, NumbOfDynVars,dyn;

  name=StringAlloc();
  rename=StringAlloc();

  NumbOfDynVars=NoOfSpec()+NumOfDynVar()-NumOfConstraint();
  for(i=1; i<=NoOfReact(); i++) {
    dyn=0;
    for(j=1; j<=NoOfSpec(); j++) {
      GetSpecNo(j, name, &charge);
      if (IsSpecInConstraint(name, charge)==0) {
        RenameSpec(rename, name, charge);
        rfwds[i-1][dyn]=TreeCreate();
        TreeDerive(rfwds[i-1][dyn], rfw[i-1], rename);
        rrvds[i-1][dyn]=TreeCreate();
        TreeDerive(rrvds[i-1][dyn], rrv[i-1], rename);
	dyn++;
      } /* if */
    } /* for j */
    for(j=0; j<NumOfDynVar(); j++) {
    /*
      NumInd=j+NoOfSpec()-NumOfConstraint();
      GetDynVarNo(j, name);
      rfwds[i][NumInd]=TreeCreate();
      TreeDerive(rfwds[i][NumInd], rfw[i], name);
      rrvds[i][NumInd]=TreeCreate();
      TreeDerive(rrvds[i][NumInd], rrv[i], name);
      */
    } /* for j */
  } /* for i */
  StringFree(name);
  StringFree(rename);
} /* GenerateDiffReac */
