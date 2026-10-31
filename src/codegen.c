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

#include <stdlib.h>

/* Global variables which can be used by the code generators. */
Tree *v;
Tree *con;
Tree **jacobi;
Tree ***hess;
Tree ****keld;
Tree *rfw;
Tree *rrv;
Tree **rfwds;
Tree **rrvds;

/***************************************************************************** 
  Code generator independent routines (exported).
  All routines are using global variables defined in the header file.
*****************************************************************************/

void InitCodeGenVar(int n, int m, int r) {

  int i, j, k;

  v=(Tree *)calloc(n, sizeof(Tree));
  con=(Tree *)calloc(m, sizeof(Tree));
  jacobi=(Tree **)calloc(n, sizeof(Tree *));
  for(i=0; i<n; i++) 
    jacobi[i]=(Tree *)calloc(n, sizeof(Tree));
  hess=(Tree ***)calloc(n, sizeof(Tree **));
  for(i=0; i<n; i++) {
    hess[i]=(Tree **)calloc(n, sizeof(Tree *));
    for(j=0; j<n; j++)
      hess[i][j]=(Tree *)calloc(n, sizeof(Tree));
  } /* for i */
  keld=(Tree ****)calloc(n, sizeof(Tree ***));
  for(i=0; i<n; i++) {
    keld[i]=(Tree ***)calloc(n, sizeof(Tree **));
    for(j=0; j<n; j++) {
      keld[i][j]=(Tree **)calloc(n, sizeof(Tree *));
      for(k=0; k<n; k++)
	keld[i][j][k]=(Tree *)calloc(n, sizeof(Tree));
    } /* for j */
  } /* for i */
  rfw=(Tree *)calloc(r, sizeof(Tree));
  rrv=(Tree *)calloc(r, sizeof(Tree));
  if(r>n) {
    rfwds=(Tree **)calloc(r, sizeof(Tree *));
    for(i=0; i<r; i++) 
      rfwds[i]=(Tree *)calloc(r, sizeof(Tree));
    rrvds=(Tree **)calloc(r, sizeof(Tree *));
    for(i=0; i<r; i++) 
      rrvds[i]=(Tree *)calloc(r, sizeof(Tree));
  } else {
    rfwds=(Tree **)calloc(n, sizeof(Tree *));
    for(i=0; i<n; i++) 
      rfwds[i]=(Tree *)calloc(n, sizeof(Tree));
    rrvds=(Tree **)calloc(n, sizeof(Tree *));
    for(i=0; i<n; i++) 
      rrvds[i]=(Tree *)calloc(n, sizeof(Tree));
  }
} /* InitCodeGenVar */


void GenerateRateExpr(void) {

  double charge, temp, coeff;
  char   *name, *rename;
  Tree   v_temp, tmp;
  Tree   *v_store, *r;
  int    i, j, react_no, finished, constraint, dyn;

  v_store=(Tree *)calloc(NoOfSpec()-NumOfConstraint()+NumOfDynVar(), sizeof(Tree));
  r=(Tree *)calloc((ExprSize+SymSize), sizeof(Tree));
  name=StringAlloc();
  rename=StringAlloc();

  /* process constraints */
  for(i=1; i<=NumOfConstraint(); i++) { 
    con[i-1]=GetConstraintNo(i, name, &charge);
    RenameSpec(rename, name, charge);
    for(j=1; j<i; j++) 
      con[j-1]=TreeSubstTree(con[j-1], rename, con[i-1]);
  } /* for i */
  for(i=1; i<=NoOfReact(); i++) {
    r[i-1]=TreeConst(0.0);
    rfw[i-1]=r[i-1];
    rrv[i-1]=r[i-1];
    react_no=GetReactNo(i-1); 
    switch (GetReactKind(react_no)) {
    case uni: 
      if (GetRateKind(react_no, uni, 1)==1) {
	r[i-1]=GetRateConst(react_no, uni, 1);
	finished=GetFirstSpecA(react_no, name, &charge, &coeff, 0);
	while (finished==1) {
	  if (coeff!=0.0) {
	    RenameSpec(rename, name, charge);
	    temp=GetPowConstInReact(react_no, name, charge, 0);
	    constraint=IsSpecInConstraint(name, charge);
	    if (constraint>0)
	      tmp=GetConstraintNo(constraint, name, &charge);
	    else
	      tmp=TreeVar(rename);
	    tmp=TreePow(tmp, TreeConst(temp));
	    r[i-1]=TreeMul(r[i-1], tmp);
	  } /* if */
	  finished=GetNextSpecA(name, &charge, &coeff, 0);
	} /* while */
      } /* if */
      else {  
	r[i-1]=GetRateExpr(react_no, uni, 1);
	for(j=1; j<=NoOfSpec(); j++) {
	  GetSpecNo(j, name, &charge);
	  constraint=IsSpecInConstraint(name, charge);
	  if (constraint>0) {
	    RenameSpec(rename, name, charge);
	    r[i-1]=TreeSubstTree(r[i-1], rename, con[constraint-1]);
	  } /* if */
	} /* for j */
      } /* else */
      rfw[i-1]=r[i-1];
      rrv[i-1]=TreeConst(0.0);
      break; /* case uni */
    case bi:  
      if (GetRateKind(react_no, bi, 1)==1) { 
	r[i-1]=GetRateConst(react_no, bi, 1);
	finished=GetFirstSpecA(react_no, name, &charge, &coeff, 0);
	while (finished==1) {
	  if (coeff!=0.0) {
	    RenameSpec(rename, name, charge);
	    temp=GetPowConstInReact(react_no, name, charge, 0);
	    constraint=IsSpecInConstraint(name, charge);
	    if (constraint>0)
	      tmp=GetConstraintNo(constraint, name, &charge);
	    else
	      tmp=TreeVar(rename);
	    tmp=TreePow(tmp, TreeConst(temp));
	    r[i-1]=TreeMul(r[i-1], tmp);
	  } /* if */
	  finished=GetNextSpecA(name, &charge, &coeff, 0);
	} /* while */
      } /* if */
      else {
	r[i-1]=GetRateExpr(react_no, bi, 1);
        for(j=1; j<=NoOfSpec(); j++) {
          GetSpecNo(j, name, &charge);
          constraint=IsSpecInConstraint(name, charge);
          if (constraint>0) {
            RenameSpec(rename, name, charge);
            r[i-1]=TreeSubstTree(r[i-1], rename, con[constraint-1]);
          } /* if */
        } /* for j */
      } /* else */
      rfw[i-1]=r[i-1];
      if (GetRateKind(react_no, bi, 2)==1) {
	v_temp=GetRateConst(react_no, bi, 2);
	finished=GetFirstSpecA(react_no, name, &charge, &coeff, 1);
	while (finished==1) {
	  if (coeff!=0.0) {
	    RenameSpec(rename, name, charge);
	    temp=GetPowConstInReact(react_no, name, charge, 1);
	    constraint=IsSpecInConstraint(name, charge);
	    if (constraint>0) 
	      tmp=GetConstraintNo(constraint, name, &charge);
	    else
	      tmp=TreeVar(rename);
	    tmp=TreePow(tmp, TreeConst(temp));
	    v_temp=TreeMul(v_temp, tmp);
	  } /* if */
	  finished=GetNextSpecA(name, &charge, &coeff, 1);
	} /* while */
      } /* if */
      else {
	v_temp=GetRateExpr(react_no, bi, 2);
        for(j=1; j<=NoOfSpec(); j++) {
          GetSpecNo(j, name, &charge);
          constraint=IsSpecInConstraint(name, charge);
          if (constraint>0) {
            RenameSpec(rename, name, charge);
            v_temp=TreeSubstTree(v_temp, rename, con[constraint-1]);
          } /* if */
        } /* for j */
      } /* else */
      r[i-1]=TreeSub(r[i-1], v_temp);
      rrv[i-1]=v_temp;
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
      tmp=TreeConst(0.0);
      for(j=1; j<=NoOfReact(); j++) {
	if (IsSpecInReact(GetReactNo(j-1), name, charge, &coeff)==1) {
	  v_temp=TreeMul(TreeConst(-coeff), r[j-1]);
	  tmp=TreeAdd(tmp, v_temp);
	} /* if */
      } /* for j */
      v_store[dyn]=tmp;
      dyn++;
    } /* if */
  } /* for i */
  for(i=1; i<=NumOfDynVar(); i++) {
    v_store[i+NoOfSpec()-NumOfConstraint()-1]=GetExprNo(i, name);
  } /* for i */
  for(i=0; i<(NoOfSpec()-NumOfConstraint()+NumOfDynVar()); i++)
    v[i]=v_store[i];

  StringFree(name);
  StringFree(rename);
} /* GenerateRateExpr */


void GenerateJacobi(void) {

  double charge;
  char   *name, *rename, *strtmp;
  int    i, j, dyn, dyn2;

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
	GetSpecNo(j, name, &charge);
	RenameSpec(rename, name, charge);
	if (IsSpecInConstraint(name, charge)==0) {
	  jacobi[dyn-1][dyn2-1]=TreeDerive(v[dyn-1], rename);
	  (void) TreeEval(jacobi[dyn-1][dyn2-1]);
	  dyn2++;
	} /* if */
      } /* for j */
      for(j=1; j<=NumOfDynVar(); j++) {
	GetDynVarNo(j, name);
	jacobi[dyn-1][j-1+NoOfSpec()-NumOfConstraint()]=TreeDerive(v[dyn-1], name);
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
	jacobi[i+NoOfSpec()-NumOfConstraint()-1][dyn-1]=TreeDerive(v[i+NoOfSpec()-NumOfConstraint()-1], rename);
      } /* if */
    } /* for j */
    for(j=1; j<=NumOfDynVar(); j++) {
      GetDynVarNo(j, name);
      jacobi[i+NoOfSpec()-NumOfConstraint()-1][j+NoOfSpec()-NumOfConstraint()-1]=TreeDerive(v[i+NoOfSpec()-NumOfConstraint()-1], name);
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
	  hess[i][j][dyn]=TreeDerive(jacobi[i][j], rename);
	  dyn++;
	} /* if */
      } /* for l */
      for(l=1; l<=NumOfDynVar(); l++) {
	GetDynVarNo(l, name);
	hess[i][j][dyn]=TreeDerive(jacobi[i][j], name);
	dyn++;
      } /* for l */
    } /* for j */
  StringFree(name);
  StringFree(rename);
} /* GenerateHessian */

void EvalJacobian(double *x, double **jac) {

  int     dyn, dyn2, i, j, l;
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
          tmp=jacobi[dyn-1][dyn2-1];
          for(l=1; l<=NoOfSpec(); l++) {
            GetSpecNo(l, name, &charge);
            RenameSpec(rename, name, charge);
            tmp=TreeSubstVar(tmp, rename, x[l-1]);
          } /* for l  */
          for(l=1; l<=NumOfDynVar(); l++) {
            GetDynVarNo(l, name);
            tmp=TreeSubstVar(tmp, name, x[l+NoOfSpec()-NumOfConstraint()-1]);
          } /* for l */
          jac[dyn-1][dyn2-1]=TreeEval(tmp);
          if (TreeGetError()==NoEval)
            fprintf(stderr, "WARNING: Was not able to evaluate element (%d, %d) of Jacobian.\n", dyn , dyn2);
          dyn2++;
        } /* if */
      } /* for j */
      dyn++;
    } /* if */
  } /* for i */
  StringFree(name);
  StringFree(rename);
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
	    keld[i][j][k][dyn]=TreeDerive(hess[i][j][k], rename);
	    dyn++;
	  } /* if */
	} /* for l */
	for(l=1; l<=NumOfDynVar(); l++) {
	  GetDynVarNo(l, name);
	  keld[i][j][k][dyn]=TreeDerive(hess[i][j][k], name);
	  dyn++;
	} /* for l */
      } /* for k */
  StringFree(name);
  StringFree(rename);
} /* GenerateKeldian */

void GenerateDiffReac(void) {

  /* It is assumed that the vectors rfw and rrv are generated. */

  double charge;
  char   *name, *rename;
  int    i, j, dyn;

  name=StringAlloc();
  rename=StringAlloc();

  for(i=1; i<=NoOfReact(); i++) {
    dyn=0;
    for(j=1; j<=NoOfSpec(); j++) {
      GetSpecNo(j, name, &charge);
      if (IsSpecInConstraint(name, charge)==0) {
        RenameSpec(rename, name, charge);
        rfwds[i-1][dyn]=TreeDerive(rfw[i-1], rename);
        rrvds[i-1][dyn]=TreeDerive(rrv[i-1], rename);
	dyn++;
      } /* if */
    } /* for j */
    for(j=0; j<NumOfDynVar(); j++) {
    /*
      NumInd=j+NoOfSpec()-NumOfConstraint();
      GetDynVarNo(j, name);
      rfwds[i][NumInd]=TreeDerive(rfw[i], name);
      rrvds[i][NumInd]=TreeDerive(rrv[i], name);
      */
    } /* for j */
  } /* for i */
  StringFree(name);
  StringFree(rename);
} /* GenerateDiffReac */
