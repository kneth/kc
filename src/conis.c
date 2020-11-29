/***************************************************************************
  KnCC - a code generator for kc and KN's continuation program.

  CopyWrong 1993-1994 by
  Kenneth Geisshirt (kneth@osc.kiku.dk)
  Department of Theoretical Chemistry
  H.C. Orsted Institute
  Universitetsparken 5
  2100 Copenhagen
  Denmark

  See kc.tex for details.

  Last updated: 11 May 1995 by KN
*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "config.h"
#include "symbmath.h"
#include "tableman.h"
#include "codegen.h"
#include "misc.h"

void ConIS(FILE *fcode) {

  double charge, temp, coeff, temp1, temp2;
  char   *name, *rename;
  time_t timer;
  Tree   v_temp, tmp, temp_tree, tmp2;
  int    i, j, k, l, react_no, finished, constraint, dyn, dyn2, dyn3, form;
  int    NumbOfParams, NumbOfDynVars;
  int    need_dd_jac, BfErrorCode;

  name=StringAlloc();
  rename=StringAlloc();
  timer=time(&timer);
  NumbOfParams=NumOfParameter();
  NumbOfDynVars=NoOfSpec()-NumOfConstraint()+NumOfDynVar();
  if ((NumbOfParams<2) || (NumbOfParams>40)) {
    fprintf(stderr, "Continuation program aborted: Wrong number of parameters - should be between 2 and 40.\n");
    return;
  } /* if */
  
  InitCodeGenVar(NoOfSpec()+NumOfDynVar()-NumOfConstraint(),
		 NumOfConstraint(), NoOfReact());
  GenerateRateExpr();
  GenerateJacobi();
  
  fprintf(fcode, "c ****************************************************\n");
  fprintf(fcode, "c  WARNING: This file was generated by kc v%s\n",
	  VERSION);
  fprintf(fcode, "c CopyWrong 1994 by Kenneth Geisshirt.\n");
  fprintf(fcode, "c %s", ctime(&timer));
  fprintf(fcode, "c ****************************************************\n");
  fprintf(fcode, "\n");
  
  /* printing derivs */
  
  fprintf(fcode, "      subroutine model(ndim,nvar,n,tvar,xvar,fvar,gvar)\n");
  fprintf(fcode, "\n");
  fprintf(fcode, "c specification of the user's problem\n");
  fprintf(fcode, "c right hand sides and jacobi matrix of\n");
  fprintf(fcode, "c the model equations are evaluated here\n");
  fprintf(fcode, "c\n");
  fprintf(fcode, "c tvar   : time (explicitly occures only for fodes)\n");
  fprintf(fcode, "c xvar() : array ndim state space variables\n");
  fprintf(fcode, "c fvar() : array ndim right hand sides depending on x;alpha,beta,par()\n");
  fprintf(fcode, "c       (in addition f depends explicitly on t for fodes)\n");
  fprintf(fcode, "c gvar(,): ndim by ndim+2 matrix of first derivatives,\n");
  fprintf(fcode, "c       gvar = [df/dx,df/dalpha,df/dbeta]\n");
  fprintf(fcode, "\n");
  fprintf(fcode, "c -------------------------------------------------------\n");
  fprintf(fcode, "      implicit real*8(a-h,o-z)\n");
  fprintf(fcode, "      dimension xvar(ndim),fvar(ndim),gvar(ndim,nvar)\n");
  fprintf(fcode, "      common/fixp/");
  for(i=3; i<=NumbOfParams; i++) {
    GetParamNo(i, name, &charge, &form);
    if (form==2)
      RenameSpec(rename, name, charge);
    else
      strcpy(rename, name);
    fprintf(fcode, "%s,", rename);
  } /* for i */
  fprintf(fcode, "dummy(%d)\n",42-NumbOfParams);
  fprintf(fcode, "      common/varp/");
  for(i=1; i<=2; i++) {
    GetParamNo(i, name, &charge, &form);
    if (form==2)
      RenameSpec(rename, name, charge);
    else
      strcpy(rename, name);
    fprintf(fcode, "%s,", rename);
  } /* for i */
  fprintf(fcode, "arg,per\n");
  fprintf(fcode, "c -------------------------------------------------------\n");
  

  dyn=1;
  for(i=1; i<=NoOfSpec(); i++) {
    GetSpecNo(i, name, &charge);
    RenameSpec(rename, name, charge);
    constraint=IsSpecInConstraint(name, charge);
    if (constraint==0) {
      fprintf(fcode, "      %s = xvar(%d)\n", rename, dyn);
      dyn++;
    } /* if */
  } /* for i*/
  fprintf(fcode, "\n");
  for(i=1; i<=NumOfDynVar(); i++) {
    GetDynVarNo(i, name);
    fprintf(fcode, "      %s = xvar(%d)\n", name, 
	    i+NoOfSpec()-NumOfConstraint());
  } /* for i */
  fprintf(fcode, "\n");
  for(i=1; i<=NoOfSpec(); i++) {
    GetSpecNo(i, name, &charge);
    RenameSpec(rename, name, charge);
    constraint=IsSpecInConstraint(name, charge);
    if (constraint>0) {
      fprintf(fcode, "      %s = ", rename);
      TreePrint(con[constraint-1], 1, fcode);
      fprintf(fcode, "\n");
    } /* if */
  } /* for i */
  fprintf(fcode, "\n");
  
  for(i=1; i<=NumbOfDynVars; i++) {
    fprintf(fcode, "      fvar(%d) = ", i);
    TreePrint(v[i-1], 1, fcode);
    fprintf(fcode, "\n");
  } /* for i */
  fprintf(fcode, "\n");
  
  for(i=1; i<=NumbOfDynVars; i++) {
    for(j=1; j<=NumbOfDynVars; j++) {
      temp=TreeEval(jacobi[i-1][j-1]);
      fprintf(fcode, "      gvar(%d,%d) = ", i, j);
      TreePrint(jacobi[i-1][j-1], 1, fcode);
      fprintf(fcode, "\n");
    } /* for j */
    fprintf(fcode, "\n");
  } /* for i */
  GetParamNo(1, name, &charge, &form);
  if (form==1)
    strcpy(rename, name);
  else
    RenameSpec(rename, name, charge);
  for(j=0; j<NumbOfDynVars; j++) {
    tmp=TreeCreate();
    TreeDerive(tmp, v[j], rename);
    fprintf(fcode, "      gvar(%d,%d) = ", j+1,NumbOfDynVars+1);
    TreePrint(tmp, 1, fcode);
    fprintf(fcode, "\n");
    TreeKill(tmp);
  } /* for j */
  fprintf(fcode, "\n\n");
  
  GetParamNo(1, name, &charge, &form);
  if (form==2)
    RenameSpec(rename, name, charge);
  else
    strcpy(rename, name);
  for(j=0; j<NumbOfDynVars; j++) {
    tmp=TreeCreate();
    TreeDerive(tmp, v[j], rename);
    fprintf(fcode, "      gvar(%d,%d) = ", j+1, NumbOfDynVars+1);
    TreePrint(tmp, 1, fcode);
    fprintf(fcode, "\n");
    TreeKill(tmp);
  } /* for j */
  fprintf(fcode, "\n");
  GetParamNo(2, name, &charge, &form);
  if (form==2)
    RenameSpec(rename, name, charge);
  else
    strcpy(rename, name);
  for(j=0; j<NumbOfDynVars; j++) {
    tmp=TreeCreate();
    TreeDerive(tmp, v[j], rename);
    fprintf(fcode, "      gvar(%d,%d) = ", j+1, NumbOfDynVars+2);
    TreePrint(tmp, 1, fcode);
    fprintf(fcode, "\n");
    TreeKill(tmp);
  } /* for i */
  fprintf(fcode, "\n\n");
  fprintf(fcode, "      return\n");
  fprintf(fcode, "      end\n\n\n");
  StringFree(name);
  StringFree(rename);
} /* ConIS */




