/**************************************************************************** 
  KinCode - a code generator for kc and KIN. 
  This is a full operating generator.   

  CopyWrong 1992-1994 by 
  Kenneth Geisshirt (kneth@osc.kiku.dk)
  Department of Theoretical Chemistry
  H.C. Orsted Institute
  Universitetsparken 5
  2100 Copenhagen 
  Denmark
  
  See kc.tex for details

  Last updated: 28 July 1994
****************************************************************************/

#include <stdio.h>
#include <malloc.h>
#include <time.h>
#include "config.h"
#include "symbmath.h"
#include "tableman.h"
#include "codegen.h"
#include "misc.h"

void KinCode(FILE *code) {

  double charge, temp, coeff;
  char   name[STRING_LENGTH], rename[STRING_LENGTH];
  time_t timer;
  Tree   v_temp, tmp, temp_tree;
  int    i, j, react_no, finished, constraint, dyn, dyn2;

  timer=time(&timer);

  InitCodeGenVar(NoOfSpec()+NumOfDynVar()-NumOfConstraint(),
		 NumOfConstraint());
  GenerateRateExpr(1, 0, 0, 0);
  GenerateJacobi(1, 0);

  fprintf(code, "(* %s*)\n", ctime(&timer));
  i=NoOfSpec()+NumOfDynVar()-NumOfConstraint(); /* abuse of i */
  fprintf(code, "CONST\n");   
  fprintf(code, "  n   = %d;\n", i);
  fprintf(code, "  np = %d;\n\n", i);
  temp=GetConstant("mode");
  if (GetError()==NoError)
    switch ((int)temp) {
    case 0:
      fprintf(code, "  PERT = FALSE;\n");
      break;
    case 1:
      fprintf(code, "  PERT = TRUE;\n");
      break;
    }
  else
    fprintf(code, "  PERT = FALSE;\n");
  fprintf(code, "VAR\n");
  for(i=1; i<=NoOfSpec(); i++) {
    GetSpecNo(i, name, &charge);
    RenameSpec(rename, name, charge);
    fprintf(code, "  %s : LONGREAL;\n", rename);
  }; /* for i */
  for(i=1; i<=NumOfDynVar(); i++) {
    GetDynVarNo(i, name);
    fprintf(code, "  %s : LONGREAL;\n", name);
  }; /* for i */
  fprintf(code, "  species : ARRAY[1..n] OF STRING[20];\n");
  fprintf(code, "\n(* FILE LIMIT *)\n\n");
  fprintf(code, "PROCEDURE derivsinit;\n");
  fprintf(code, "BEGIN\n");
  dyn=1;
  for(i=1; i<=NoOfSpec(); i++) {
    GetSpecNo(i, name, &charge);
    constraint=IsSpecInConstraint(name, charge);
    if (constraint==0) {
      if (charge==0.0)
        fprintf(code, "  species[%d]:='%s';\n", dyn, name);
      else { 
	if (charge==MAXFLOAT)
	  fprintf(code, "  species[%d]:='%s.';\n", dyn, name);
	else { 
	  if (charge>0.0) {
	  fprintf(code, "  species[%d]:='%s", dyn, name);
	  for(j=1; j<=(int)charge; j++)
	    fprintf(code, "+");              
          fprintf(code, "';\n");
        } /* if */ 
	  else {         
	    fprintf(code, "  species[%d]:='%s", dyn, name);
	    for(j=1;j<=-(int)charge;j++)
	      fprintf(code, "-");
	    fprintf(code, "';\n");
	  }; /* else */
	}; /* else */
      }; /* else */
      dyn++;
    }; /* if */
  }; /* for i */
  for(i=1; i<=NumOfDynVar(); i++) {
    GetDynVarNo(i, name);
    fprintf(code, "  species[%d]:='%s';\n", i+NoOfSpec(), name);
  }; /* for i */

  GetAndPrintConst("stime", "tb", 1, 1.0, code, 2);
  GetAndPrintConst("dtime", "dt", 1, 1.0, code, 2);
  GetAndPrintConst("etime", "te", 1, 10.0, code, 2);
  GetAndPrintConst("htime", "hb", 1, 1.0, code, 2);
  GetAndPrintConst("epsr", "epsr", 1, 1e-3, code, 2);
  GetAndPrintConst("epsa", "epsa", 1, 1e-20, code, 2);
  GetAndPrintConst("ptime", "perttime", 1, 1.0, code, 2);
  GetStrConst("name", name);
  if (GetError()==NotFound)
    fprintf(code, "name_datafile:='kinwrk.dat';\n\n");
  else
    fprintf(code, "name_datafile:='%s';\n\n", name);
  dyn=1; 
  for(i=1; i<=NoOfSpec(); i++) {
    GetSpecNo(i, name, &charge);
    if (IsSpecInConstraint(name, charge)==0) { 
      temp=GetBeginConc(name, charge);
      fprintf(code, "  xx[%d]:=%e;\n", dyn, temp);
      dyn++;
    }; /* if */
  }; /* for i */
  dyn=1;
  for(i=1; i<=NoOfSpec(); i++) {
    GetSpecNo(i, name, &charge);
    if (IsSpecInConstraint(name, charge)==0) {
      temp=GetSpecConst(name, charge, "pert");
      if (GetError()==NoError)
	fprintf(code, "xxpert[%d]:=%e;\n", dyn, temp);
      else
	fprintf(code, "xxpert[%d]:=0.0;\n", dyn);
      dyn++;
    } /* if */
  } /* for i */
  for(i=1; i<=NumOfExpr(); i++) {
    GetDynVarNo(i, name);
    fprintf(code, "  xx[%d]:=%e;\n", i+NoOfSpec()-NumOfConstraint(), GetInitValue(name));
  } /* for i */
  fprintf(code, "END;\n\n");
  fprintf(code, "PROCEDURE derivs(bj:BOOLEAN; xx_:glnarray; t:LONGREAL;\n");
  fprintf(code, "               VAR vv_:glnarray; VAR jj_:glnpbynp);\n");
  fprintf(code, "BEGIN\n");
  dyn=1;
  for(i=1; i<=NoOfSpec(); i++) {
    GetSpecNo(i, name, &charge);
    RenameSpec(rename, name, charge);
    constraint=IsSpecInConstraint(name, charge);
    if (constraint==0) {  
      fprintf(code, "%s:=xx_[%d];\n", rename, dyn);
      dyn++;
    }; /* if */
  }; /* for i*/
  for(i=1; i<=NumOfDynVar(); i++) {
    GetDynVarNo(i, name);
    fprintf(code, "%s:=xx_[%d];\n", name, i+NoOfSpec()-NumOfConstraint());
  }; /* for i */
  for(i=1; i<=NoOfSpec(); i++) {
    GetSpecNo(i, name, &charge);
    RenameSpec(rename, name, charge);
    constraint=IsSpecInConstraint(name, charge);
    if (constraint>0) {  
      fprintf(code, "%s:=", rename);
      TreePrint(con[constraint-1], 2, code);
      fprintf(code, ";\n");
    }; /* if */
  }; /* for i */

  for(i=1; i<=(NoOfSpec()-NumOfConstraint()+NumOfDynVar()); i++) {
    fprintf(code, "vv_[%d]:=", i);
    TreePrint(v[i-1], 2, code);
    fprintf(code, ";\n");
  } /* for i */

  fprintf(code, "IF bj THEN BEGIN\n");
  for(i=1; i<=(NoOfSpec()-NumOfConstraint()+NumOfDynVar()); i++) 
    for(j=1; j<=(NoOfSpec()-NumOfConstraint()+NumOfDynVar()); j++) {
      temp=TreeEval(jacobi[i-1][j-1]);
      if (TreeGetError()==NoEval) {
	fprintf(code, "jj_[%d,%d]:=", i, j);
	TreePrint(jacobi[i-1][j-1], 2, code);
	fprintf(code, ";\n");
      } /* if */
      else
	if (temp!=0.0)
	  fprintf(code, "jj_[%d,%d]:=%e;\n", i, j, temp);
    } /* for j */
  fprintf(code, "  END\nEND;\n");
} /* KinCode */
