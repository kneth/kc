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
#include <malloc.h>
#include <time.h>
#include <math.h>
#include "config.h"
#include "symbmath.h"
#include "tableman.h"
#include "codegen.h"
#include "misc.h"

void KnCC(FILE *ccode, FILE *hcode) {

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
  if ((NumbOfParams==0) || (NumbOfParams>2)) {
    fprintf(stderr, "KNCont: Wrong number of parameters - should be either 1 or 2.\n");
    return;
  } /* if */

  InitCodeGenVar(NoOfSpec()+NumOfDynVar()-NumOfConstraint(),
		 NumOfConstraint(), NoOfReact());
  GenerateRateExpr();
  GenerateJacobi();

  fprintf(ccode, "c ****************************************************\n");
  fprintf(ccode, "c  WARNING: This file was generated by kc v%s\n",
	  VERSION);
  fprintf(ccode, "c CopyWrong 1994 by Kenneth Geisshirt.\n");
  fprintf(ccode, "c %s", ctime(&timer));
  fprintf(ccode, "c ****************************************************\n");
  fprintf(ccode, "\n");

  /* printing derivs */

  fprintf(ccode, "      subroutine model(ndim,nvar,n,t,x,f,g)\n");
  fprintf(ccode, "\n");
  fprintf(ccode, "c specification of the user's problem\n");
  fprintf(ccode, "c right hand sides and jacobi matrix of\n");
  fprintf(ccode, "c the model equations are evaluated here\n");
  fprintf(ccode, "c\n");
  fprintf(ccode, "c t   : time (explicitly occures only for fodes)\n");
  fprintf(ccode, "c x() : array ndim state space variables\n");
  fprintf(ccode, "c f() : array ndim right hand sides depending on x;alpha,beta,par()\n");
  fprintf(ccode, "c       (in addition f depends explicitly on t for fodes)\n");
  fprintf(ccode, "c g(,): ndim by ndim+2 matrix of first derivatives,\n");
  fprintf(ccode, "c       g = [df/dx,df/dalpha,df/dbeta]\n");
  fprintf(ccode, "\n");
  fprintf(ccode, "c -------------------------------------------------------\n");
  fprintf(ccode, "      implicit real*8(a-h,o-z)\n");
  fprintf(ccode, "      dimension x(ndim),f(ndim),g(ndim,nvar)\n");
  fprintf(ccode, "      common/fixp/dummy(40)\n");
  fprintf(ccode, "      common/varp/");
  for(i=1; i<=NumbOfParams; i++) {
    GetParamNo(i, name, &charge, &form);
    if (form==2)
      RenameSpec(rename, name, charge);
    else
      strcpy(rename, name);
    fprintf(ccode, "%s,", rename);
  } /* for i */
  fprintf(ccode, "arg,per\n");
  fprintf(ccode, "c -------------------------------------------------------\n");


  for(i=1; i<=NoOfSpec(); i++) {
    GetSpecNo(i, name, &charge);
    RenameSpec(rename, name, charge);
    fprintf(ccode, "      double %s\n", rename);
  } /* for i */
  for(i=1; i<=NumOfDynVar(); i++) {
    GetDynVarNo(i, name);
    fprintf(ccode, "double %s\n", name);
  } /* for i */

  dyn=1;
  for(i=1; i<=NoOfSpec(); i++) {
    GetSpecNo(i, name, &charge);
    RenameSpec(rename, name, charge);
    constraint=IsSpecInConstraint(name, charge);
    if (constraint==0) {
      fprintf(ccode, "      %s = x(%d)\n", rename, dyn);
      dyn++;
    } /* if */
  } /* for i*/
  fprintf(ccode, "\n");
  for(i=1; i<=NumOfDynVar(); i++) {
    GetDynVarNo(i, name);
    fprintf(ccode, "      %s = x(%d)\n", name, i+NoOfSpec()-NumOfConstraint());
  } /* for i */
  fprintf(ccode, "\n");
  for(i=1; i<=NoOfSpec(); i++) {
    GetSpecNo(i, name, &charge);
    RenameSpec(rename, name, charge);
    constraint=IsSpecInConstraint(name, charge);
    if (constraint>0) {
      fprintf(ccode, "      %s = ", rename);
      TreePrint(con[constraint-1], 1, ccode);
      fprintf(ccode, "\n");
    }; /* if */
  }; /* for i */
  fprintf(ccode, "\n");

  for(i=1; i<=NumbOfDynVars; i++) {
    fprintf(ccode, "      f(%d) = ", i);
    TreePrint(v[i-1], 1, ccode);
    fprintf(ccode, "\n");
  } /* for i */
  fprintf(ccode, "\n");

  for(i=1; i<=NumbOfDynVars; i++) {
    for(j=1; j<=NumbOfDynVars; j++) {
      temp=TreeEval(jacobi[i-1][j-1]);
  	fprintf(ccode, "      g(%d,%d) = ", i, j);
	TreePrint(jacobi[i-1][j-1], 1, ccode);
	fprintf(ccode, "\n");
    } /* for j */
    fprintf(ccode, "\n");
  } /* for i */
    GetParamNo(1, name, &charge, &form);
    if (form==1)
      strcpy(rename, name);
    else
      RenameSpec(rename, name, charge);
    for(j=0; j<NumbOfDynVars; j++) {
      tmp=TreeCreate();
      TreeDerive(tmp, v[j], rename);
      fprintf(ccode, "      g(%d,%d) = ", j+1,NumbOfDynVars+1);
      TreePrint(tmp, 1, ccode);
      fprintf(ccode, "\n");
      TreeKill(tmp);
    } /* for j */
    fprintf(ccode, "\n\n");

    GetParamNo(1, name, &charge, &form);
    if (form==2)
      RenameSpec(rename, name, charge);
    else
      strcpy(rename, name);
    for(j=0; j<NumbOfDynVars; j++) {
      tmp=TreeCreate();
      TreeDerive(tmp, v[j], rename);
      fprintf(ccode, "      g(%d,%d) = ", j+1, NumbOfDynVars+1);
      TreePrint(tmp, 1, ccode);
      fprintf(ccode, "\n");
      TreeKill(tmp);
    } /* for j */
    fprintf(ccode, "\n");
    GetParamNo(2, name, &charge, &form);
    if (form==2)
      RenameSpec(rename, name, charge);
    else
      strcpy(rename, name);
    for(j=0; j<NumbOfDynVars; j++) {
      tmp=TreeCreate();
      TreeDerive(tmp, v[j], rename);
      fprintf(ccode, "      g(%d,%d) = ", j+1, NumbOfDynVars+2);
      TreePrint(tmp, 1, ccode);
      fprintf(ccode, "\n");
      TreeKill(tmp);
    } /* for i */
    fprintf(ccode, "\n\n");
  fprintf(ccode, "      return\n");
  fprintf(ccode, "      end\n\n\n");
  StringFree(name);
  StringFree(rename);
} /* KnCC */
