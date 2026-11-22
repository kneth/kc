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
expr_t *v;
expr_t *con;
expr_t **jacobi;
expr_t ***hess;
expr_t ****keld;
expr_t *rfw;
expr_t *rrv;
expr_t **rfwds;
expr_t **rrvds;

/***************************************************************************** 
  Code generator independent routines (exported).
  All routines are using global variables defined in the header file.
*****************************************************************************/

void InitCodeGenVar(int n, int m, int r) {

  int i, j, k;

  v=(expr_t *)calloc(n, sizeof(expr_t));
  con=(expr_t *)calloc(m, sizeof(expr_t));
  jacobi=(expr_t **)calloc(n, sizeof(expr_t *));
  for(i=0; i<n; i++) 
    jacobi[i]=(expr_t *)calloc(n, sizeof(expr_t));
  hess=(expr_t ***)calloc(n, sizeof(expr_t **));
  for(i=0; i<n; i++) {
    hess[i]=(expr_t **)calloc(n, sizeof(expr_t *));
    for(j=0; j<n; j++)
      hess[i][j]=(expr_t *)calloc(n, sizeof(expr_t));
  } /* for i */
  keld=(expr_t ****)calloc(n, sizeof(expr_t ***));
  for(i=0; i<n; i++) {
    keld[i]=(expr_t ***)calloc(n, sizeof(expr_t **));
    for(j=0; j<n; j++) {
      keld[i][j]=(expr_t **)calloc(n, sizeof(expr_t *));
      for(k=0; k<n; k++)
	keld[i][j][k]=(expr_t *)calloc(n, sizeof(expr_t));
    } /* for j */
  } /* for i */
  rfw=(expr_t *)calloc(r, sizeof(expr_t));
  rrv=(expr_t *)calloc(r, sizeof(expr_t));
  if(r>n) {
    rfwds=(expr_t **)calloc(r, sizeof(expr_t *));
    for(i=0; i<r; i++) 
      rfwds[i]=(expr_t *)calloc(r, sizeof(expr_t));
    rrvds=(expr_t **)calloc(r, sizeof(expr_t *));
    for(i=0; i<r; i++) 
      rrvds[i]=(expr_t *)calloc(r, sizeof(expr_t));
  } else {
    rfwds=(expr_t **)calloc(n, sizeof(expr_t *));
    for(i=0; i<n; i++) 
      rfwds[i]=(expr_t *)calloc(n, sizeof(expr_t));
    rrvds=(expr_t **)calloc(n, sizeof(expr_t *));
    for(i=0; i<n; i++) 
      rrvds[i]=(expr_t *)calloc(n, sizeof(expr_t));
  }
} /* InitCodeGenVar */


void GenerateRateExpr(void) {

  double charge, temp, coeff;
  char   *name, *rename;
  expr_t v_temp, tmp;
  expr_t *v_store, *r;
  int    i, j, react_no, finished, constraint, dyn;

  v_store=(expr_t *)calloc(NoOfSpec()-NumOfConstraint()+NumOfDynVar(), sizeof(expr_t));
  r=(expr_t *)calloc((ExprSize+SymSize), sizeof(expr_t));
  name=StringAlloc();
  rename=StringAlloc();

  /* process constraints */
  for(i=1; i<=NumOfConstraint(); i++) { 
    con[i-1]=GetConstraintNo(i, name, &charge);
    RenameSpec(rename, name, charge);
    for(j=1; j<i; j++) 
      con[j-1]=expr_subst_expr(con[j-1], rename, con[i-1]);
  } /* for i */
  for(i=1; i<=NoOfReact(); i++) {
    r[i-1]=expr_const(0.0);
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
	      tmp=expr_var(rename);
	    tmp=expr_pow(tmp, expr_const(temp));
	    r[i-1]=expr_mul(r[i-1], tmp);
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
	    r[i-1]=expr_subst_expr(r[i-1], rename, con[constraint-1]);
	  } /* if */
	} /* for j */
      } /* else */
      rfw[i-1]=r[i-1];
      rrv[i-1]=expr_const(0.0);
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
	      tmp=expr_var(rename);
	    tmp=expr_pow(tmp, expr_const(temp));
	    r[i-1]=expr_mul(r[i-1], tmp);
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
            r[i-1]=expr_subst_expr(r[i-1], rename, con[constraint-1]);
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
	      tmp=expr_var(rename);
	    tmp=expr_pow(tmp, expr_const(temp));
	    v_temp=expr_mul(v_temp, tmp);
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
            v_temp=expr_subst_expr(v_temp, rename, con[constraint-1]);
          } /* if */
        } /* for j */
      } /* else */
      r[i-1]=expr_sub(r[i-1], v_temp);
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
      tmp=expr_const(0.0);
      for(j=1; j<=NoOfReact(); j++) {
	if (IsSpecInReact(GetReactNo(j-1), name, charge, &coeff)==1) {
	  v_temp=expr_mul(expr_const(-coeff), r[j-1]);
	  tmp=expr_add(tmp, v_temp);
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

  free(v_store);
  free(r);
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
	  jacobi[dyn-1][dyn2-1]=expr_derive(v[dyn-1], rename);
	  (void) expr_eval(jacobi[dyn-1][dyn2-1]);
	  dyn2++;
	} /* if */
      } /* for j */
      for(j=1; j<=NumOfDynVar(); j++) {
	GetDynVarNo(j, name);
	jacobi[dyn-1][j-1+NoOfSpec()-NumOfConstraint()]=expr_derive(v[dyn-1], name);
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
	jacobi[i+NoOfSpec()-NumOfConstraint()-1][dyn-1]=expr_derive(v[i+NoOfSpec()-NumOfConstraint()-1], rename);
      } /* if */
    } /* for j */
    for(j=1; j<=NumOfDynVar(); j++) {
      GetDynVarNo(j, name);
      jacobi[i+NoOfSpec()-NumOfConstraint()-1][j+NoOfSpec()-NumOfConstraint()-1]=expr_derive(v[i+NoOfSpec()-NumOfConstraint()-1], name);
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
	  hess[i][j][dyn]=expr_derive(jacobi[i][j], rename);
	  dyn++;
	} /* if */
      } /* for l */
      for(l=1; l<=NumOfDynVar(); l++) {
	GetDynVarNo(l, name);
	hess[i][j][dyn]=expr_derive(jacobi[i][j], name);
	dyn++;
      } /* for l */
    } /* for j */
  StringFree(name);
  StringFree(rename);
} /* GenerateHessian */

void EvalJacobian(double *x, double **jac) {

  int     dyn, dyn2, i, j, l;
  expr_t  tmp;
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
            tmp=expr_subst_var(tmp, rename, x[l-1]);
          } /* for l  */
          for(l=1; l<=NumOfDynVar(); l++) {
            GetDynVarNo(l, name);
            tmp=expr_subst_var(tmp, name, x[l+NoOfSpec()-NumOfConstraint()-1]);
          } /* for l */
          jac[dyn-1][dyn2-1]=expr_eval(tmp);
          if (expr_get_error()==EXPR_NO_EVAL)
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
	    keld[i][j][k][dyn]=expr_derive(hess[i][j][k], rename);
	    dyn++;
	  } /* if */
	} /* for l */
	for(l=1; l<=NumOfDynVar(); l++) {
	  GetDynVarNo(l, name);
	  keld[i][j][k][dyn]=expr_derive(hess[i][j][k], name);
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
        rfwds[i-1][dyn]=expr_derive(rfw[i-1], rename);
        rrvds[i-1][dyn]=expr_derive(rrv[i-1], rename);
	dyn++;
      } /* if */
    } /* for j */
    for(j=0; j<NumOfDynVar(); j++) {
    /*
      NumInd=j+NoOfSpec()-NumOfConstraint();
      GetDynVarNo(j, name);
      rfwds[i][NumInd]=expr_derive(rfw[i], name);
      rrvds[i][NumInd]=expr_derive(rrv[i], name);
      */
    } /* for j */
  } /* for i */
  StringFree(name);
  StringFree(rename);
} /* GenerateDiffReac */
