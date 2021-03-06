/************************************************************************
  A code generator which is using VODE [1] as integrator.

  CopyWrong 1994 by
  Kenneth Geisshirt (kneth@osc.kiku.dk)
  Department of Theoretical Chemistry
  H.C. Orsted Institute
  Universitetsparken 5
  2100 Copenhagen
  Denmark

  Reference:
  [1] VODE: A Variable-Coefficient ODE Solver.
      P.N. Brown et al., SIAM J. Sci. Stat. Comp.
      pp. 1038-1051, volume 10, number 5, 1989.

  Last updated: 17 August 1994
************************************************************************/

#include <stdio.h>
#include <malloc.h>
#include <time.h>

#include "config.h"
#include "tableman.h"
#include "symbmath.h"
#include "misc.h"
#include "codegen.h"

void VODE(FILE *code) {

  double    charge, temp, coeff;
  char      name[STRING_LENGTH], rename[STRING_LENGTH];
  time_t    timer;
  Tree      v_tree, tmp, temp_tree;
  int       i, j, react_no, finished, constraint, dyn, dyn2, NumbDynVars;

  timer=time(&timer);

  InitCodeGenVar(NoOfSpec()-NumOfConstraint()+NumOfDynVar(), 
		 NumOfConstraint());
  
  fprintf(code, "C--------------------------------------------------------\n");
  fprintf(code, "C WARNING: This file was generated by kc v%s\n", VERSION);
  fprintf(code, "C CopyWrong Kenneth Geisshirt\n");
  fprintf(code, "C %s", ctime(&timer));
  fprintf(code, "C--------------------------------------------------------\n");

  NumbDynVars=NoOfSpec()+NumOfDynVar()-NumOfConstraint();
  fprintf(code, "      PROGRAM ODESOLV\n");
  fprintf(code, "      DOUBLE PRECISION ATOL,RPAR,RTOL,RWORK,T,TOUT,Y\n");
  fprintf(code, "      DIMENSION Y(%d),RWORK(%d),IWORK(%d),RPAR(%d),IPAR(%d)\n", 
	  NumbDynVars, 22+9*NumbDynVars+2*NumbDynVars*NumbDynVars,
	  30+NumbDynVars, 22+9*NumbDynVars+2*NumbDynVars*NumbDynVars, 
	  30+NumbDynVars);
  fprintf(code, "      INTEGER ITOL,ITASK,ISTATE,IOPT,LRW,LIW,MF,I,NEQ\n");
  fprintf(code, "      NEQ=%d\n", NumbDynVars);
  dyn=1;
  for(i=1; i<=NoOfSpec(); i++) {
    GetSpecNo(i, name, &charge);
    if (IsSpecInConstraint(name, charge)==0) {
      temp=GetBeginConc(name, charge);
      if (GetError()==NotFound)
	fprintf(code, "      Y(%d)=0.0\n", dyn);
      else
	fprintf(code, "      Y(%d)=%e\n", dyn, temp);
      dyn++;
    }
  }
  for(i=1; i<=NumOfDynVar(); i++) {
    GetDynVarNo(i, name);
    temp=GetInitValue(name);
    if (GetError()==NotFound)
      fprintf(code, "      Y(%d)=0.0\n", i+NoOfSpec()-NumOfConstraint());
    else
      fprintf(code, "      Y(%d)=%e\n", i+NoOfSpec()-NumOfConstraint(), temp);
  }
  fprintf(code, "      T=0.0D0\n");
  fprintf(code, "      ITOL=1\n");
  fprintf(code, "      ITASK=1\n");
  fprintf(code, "      ISTATE=1\n");
  fprintf(code, "      IOPT=0\n");
  fprintf(code, "      LRW=%d\n", 22+9*NumbDynVars+2*NumbDynVars*NumbDynVars);
  fprintf(code, "      LIW=%d\n", 30+NumbDynVars);
  fprintf(code, "      MF=21\n");
  GetAndPrintConst("epsr", "RTOL", 1, 1e-10, code, 1);
  GetAndPrintConst("epsa", "ATOL", 1, 1e-10, code, 1);
  GetAndPrintConst("dtime", "TOUT", 1, 1.0, code, 1);
  fprintf(code, "      PRINT *,'BEGINNING INTEGRATION'\n");
  fprintf(code, " 1000 CALL DVODE(F,NEQ,Y,T,TOUT,ITOL,RTOL,ATOL,ITASK,\n");
  fprintf(code, "     &     ISTATE,IOPT,RWORK,LRW,IWORK,LIW,JAC,\n");
  fprintf(code, "     &     MF,RPAR,IPAR)\n");
  fprintf(code, "      PRINT *,T,(Y(I),I=1,NEQ)\n");
  temp=GetConstant("etime");
  if (GetError()==NotFound)
    temp=10.0;
  fprintf(code, "      IF (T.GT.%e) GOTO 1001\n", temp);
  temp=GetConstant("dtime");
  if (GetError()==NotFound)
    temp=1.0;
  fprintf(code, "      TOUT=TOUT+%e\n", temp);
  fprintf(code, "      GOTO 1000\n");
  fprintf(code, " 1001 END\n");

  GenerateRateExpr(1, 0, 0, 0);
  GenerateJacobi(1, 0);

  fprintf(code, "      SUBROUTINE F(NEQ,T,YSTATE,YDOT,RPAR,IPAR)\n");
  fprintf(code, "      DOUBLE PRECISION T,YSTATE,YDOT,RPAR\n");
  fprintf(code, "      DIMENSION YSTATE(NEQ),YDOT(NEQ)\n");
  for(i=1; i<=NoOfSpec(); i++) {
    GetSpecNo(i, name, &charge);
    if (IsSpecInConstraint(name, charge)==0) {
      RenameSpec(rename, name, charge);
      fprintf(code, "      DOUBLE PRECISION %s\n", rename);
    }
  }
  for(i=1; i<=NumOfDynVar(); i++) {
    GetDynVarNo(i, name);
    fprintf(code, "      DOUBLE PRECISION %s\n", name);
  }
  dyn=1;
  for(i=1; i<=NoOfSpec(); i++) {
    GetSpecNo(i, name, &charge);
    if (IsSpecInConstraint(name, charge)==0) {
      RenameSpec(rename, name, charge);
      fprintf(code, "      %s=YSTATE(%d)\n", rename, dyn);
      dyn++;
    } /* if */
  } /* for i */
  for(i=1; i<=NumOfDynVar(); i++) {
    GetDynVarNo(i, name);
    fprintf(code, "      %s=YSTATE(%d)\n", name, i+NoOfSpec()-NumOfConstraint());
  } /* for i */
  for(i=0; i<NumbDynVars; i++) {
    fprintf(code, "      YDOT(%d)=", i+1);
    TreePrint(v[i], 1, code);
    fprintf(code, "\n");
  }
  fprintf(code, "      END\n");
  fprintf(code, "      SUBROUTINE JAC(NEQ,T,YSTATE,ML,MU,PD,NROWPD,RPAR,IPAR)\n");
  fprintf(code, "      DOUBLE PRECISION T,YSTATE,PD,RPAR\n");
  fprintf(code, "      DIMENSION YSTATE(NEQ),PD(NROWPD,NEQ)\n");
  for(i=1; i<=NoOfSpec(); i++) {
    GetSpecNo(i, name, &charge);
    if (IsSpecInConstraint(name, charge)==0) {
      RenameSpec(rename, name, charge);
      fprintf(code, "      DOUBLE PRECISION %s\n", rename);
    }
  }
  for(i=1; i<=NumOfDynVar(); i++) {
    GetDynVarNo(i, name);
    fprintf(code, "      DOUBLE PRECISION %s\n", name);
  }
  dyn=1;
  for(i=1; i<=NoOfSpec(); i++) {
    GetSpecNo(i, name, &charge);
    if (IsSpecInConstraint(name, charge)==0) {
      RenameSpec(rename, name, charge);
      fprintf(code, "      %s=YSTATE(%d)\n", rename, dyn);
      dyn++;
    } /* if */
  } /* for i */
  for(i=1; i<=NumOfDynVar(); i++) {
    GetDynVarNo(i, name);
    fprintf(code, "      %s=YSTATE(%d)\n", name, i+NoOfSpec()-NumOfConstraint());
  } /* for i */
  for(i=1; i<=NumbDynVars; i++) {
    for(j=1; j<=NumbDynVars; j++) {
      fprintf(code, "      PD(%d,%d)=", i, j);
      TreePrint(jacobi[i-1][j-1], 1, code);
      fprintf(code, "\n");
    }
  }
  fprintf(code, "      END\n");
} /* VODE */



