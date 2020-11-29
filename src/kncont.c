/***************************************************************************
  KnCont - a code generator for kc and KN's continuation program.

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

void KnCont(FILE *pcode, FILE *hcode) {

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
  if (NumbOfParams==2) {
    GenerateHessian();
  } /* if */

  fprintf(hcode, "(******************************************************\n");
  fprintf(hcode, "  WARNING: This file was generated by kc v%s\n",
	  VERSION);
  fprintf(hcode, "  CopyWrong 1994 by Kenneth Geisshirt.\n");
  fprintf(hcode, "  %s", ctime(&timer));
  fprintf(hcode, "*******************************************************)\n");
  fprintf(hcode, "CONST\n");
  fprintf(hcode, "n=%d;\n", NumbOfDynVars);
  fprintf(hcode, "np=%d;\n", NumbOfDynVars);
  fprintf(hcode, "VAR\n");
  for(i=1; i<=NoOfSpec(); i++) {
    GetSpecNo(i, name, &charge);
    RenameSpec(rename, name, charge);
    fprintf(hcode, "%s:LONGREAL;\n", rename);
  } /* for i */
  for(i=1; i<=NumOfDynVar(); i++) {
    GetDynVarNo(i, name);
    fprintf(hcode, "%s:LONGREAL;\n", name);
  } /* for i */
  for(i=1; i<=NumOfParameter(); i++) {
    GetParamNo(i, name, &charge, &form);
    if (form==1)
      strcpy(rename, name);
    else
      RenameSpec(rename, name, charge);
    fprintf(hcode, "%s:LONGREAL;\n", rename);
  } /* for i */
  fprintf(hcode, "species : ARRAY[1..n] OF STRING[20];\n");

  fprintf(pcode, "(******************************************************\n");
  fprintf(pcode, "  WARNING: This file was generated by kc v%s\n",
	  VERSION);
  fprintf(pcode, "  CopyWrong 1994 by Kenneth Geisshirt.\n");
  fprintf(pcode, "  %s", ctime(&timer));
  fprintf(pcode, "*******************************************************)\n");

  /* printing Hessian */
  fprintf(pcode, "PROCEDURE djacobian(xx_:glnarray; VAR dS_:gldjacobian);\n");
  fprintf(pcode, "BEGIN\n");
  switch (NumbOfParams) {
  case 1:
    need_dd_jac=0;
    break;
  case 2:
    dyn=1;
    for(i=1; i<=NoOfSpec(); i++) {
      GetSpecNo(i, name, &charge);
      if ((IsSpecInConstraint(name, charge)==0) &&
	  (IsSpecParam(name, charge)!=1)) {
	RenameSpec(rename, name, charge);
	fprintf(pcode, "%s:=xx_[%d];\n", rename, dyn);
	dyn++;
      } /* if */
    } /* for i */
    for(i=1; i<=NumOfDynVar(); i++) {
      GetDynVarNo(i, name);
      fprintf(pcode, "%s:=xx_[%d];\n", rename,
	      i+NoOfSpec()/NumOfConstraint());
    } /* for i */
    need_dd_jac=0;

/*
    for(i=0; i<NumbOfDynVars; i++) {
      fprintf(pcode, "dS_[%d,%d,%d]:=", i+1, i+1, i+1);
      temp=TreeEval(hess[i][i][i]);
      if (TreeGetError()==NoEval) {
	need_dd_jac=1;
	TreePrint(hess[i][i][i], 2, pcode);
      } else
	fprintf(pcode, "%e", temp);
      fprintf(pcode, ";\n");
    }*/ /* for i */

    for(i=0; i<NumbOfDynVars; i++)
      for(j=0; j<NumbOfDynVars; j++)
      /*
	for(l=0; l<=j; l++) {
	*/
	for(l=0; l<NumbOfDynVars; l++) {
	  fprintf(pcode, "dS_[%d,%d,%d]:=", i+1, j+1, l+1);
	  temp=TreeEval(hess[i][j][l]);
	  if (TreeGetError()==NoEval) {
	    need_dd_jac=1;
	    TreePrint(hess[i][j][l], 2, pcode);
	  } else
	    fprintf(pcode, "%e", temp);
	  fprintf(pcode, ";\n");
	  /*
	  if (l<j) {
	  fprintf(pcode, "dS_[%d,%d,%d]:=dS_[%d,%d,%d];\n", i+1, l+1,
		  j+1, i+1, j+1, l+1);
          }*/ /* if */
	} /* for l */
    break;
  } /* switch NumbOfParams */
  fprintf(pcode, "END; (* djacobian *)\n\n");

  /* End of printing the Hessian */

  /* printing Keldian */
  fprintf(pcode, "PROCEDURE ddjacobian(xx_:glnarray;VAR ddS_:glddjacobian);\n");
  fprintf(pcode, "BEGIN\n");
  if (need_dd_jac==1) {

  if (NumbOfParams==2) {
    GenerateKeldian();
  } /* if */

  dyn=1;
  for(i=1; i<=NoOfSpec(); i++) {
    GetSpecNo(i, name, &charge);
    if ((IsSpecInConstraint(name, charge)==0) &&
	(IsSpecParam(name, charge)!=1)) {
      RenameSpec(rename, name, charge);
      fprintf(pcode, "%s:=xx_[%d];\n", rename, dyn);
      dyn++;
    } /* if */
  } /* for i */
  for(i=1; i<=NumOfDynVar(); i++) {
    GetDynVarNo(i, name);
    fprintf(pcode, "%s:=xx_[%d];\n", name, i+NoOfSpec()-NumOfConstraint());
  } /* for i */
    /*
    for(i=0; i<NumbOfDynVars; i++) {
      fprintf(pcode, "dds_[%d,%d,%d,%d]:=", i+1, i+1, i+1, i+1);
      TreePrint(keld[i][i][i][i], 2, pcode);
      fprintf(pcode, ";\n");
    }*/ /* i */


    for(i=0; i<NumbOfDynVars; i++)
      for(j=0; j<NumbOfDynVars; j++)
	for(l=0; l<NumbOfDynVars; l++)
	/*
	  for(k=0; k<=l; k++) {
	  */
	  for(k=0; k<NumbOfDynVars; k++) {
	    fprintf(pcode, "dds_[%d,%d,%d,%d]:=", i+1, j+1, l+1, k+1);
	    TreePrint(keld[i][j][l][k], 2, pcode);
	    fprintf(pcode, ";\n");
	    /*
	    fprintf(pcode, "dds_[%d,%d,%d,%d]:=ddS_[%d,%d,%d,%d];\n",
		    i+1, j+1, k+1, l+1, i+1, j+1, l+1, k+1);
	    fprintf(pcode, "dds_[%d,%d,%d,%d]:=ddS_[%d,%d,%d,%d];\n",
		    i+1, l+1, j+1, k+1, i+1, j+1, l+1, k+1);
	    fprintf(pcode, "dds_[%d,%d,%d,%d]:=ddS_[%d,%d,%d,%d];\n",
		    i+1, l+1, k+1, j+1, i+1, j+1, l+1, k+1);
	    fprintf(pcode, "dds_[%d,%d,%d,%d]:=ddS_[%d,%d,%d,%d];\n",
		    i+1, k+1, l+1, j+1, i+1, j+1, l+1, k+1);
	    fprintf(pcode, "dds_[%d,%d,%d,%d]:=ddS_[%d,%d,%d,%d];\n",
		    i+1, k+1, j+1, l+1, i+1, j+1, l+1, k+1);
	  */
	  } /* for k */
  } /* need_dd_jac==1 */
  fprintf(pcode, "END; (* ddjacobian *)\n\n");

  /* End of printing the Keldian */


  /* printing derivs */

  fprintf(pcode, "PROCEDURE derivs(bj_:BOOLEAN; xx_:glnarray; t_:LONGREAL;\n");
  fprintf(pcode, "                 VAR vv_:glnarray; VAR jj_:glnpbynp);\n");
  fprintf(pcode, "BEGIN\n");
  dyn=1;
  for(i=1; i<=NoOfSpec(); i++) {
    GetSpecNo(i, name, &charge);
    RenameSpec(rename, name, charge);
    constraint=IsSpecInConstraint(name, charge);
    if (constraint==0) {
      fprintf(pcode, "%s:=xx_[%d];\n", rename, dyn);
      dyn++;
    } /* if */
  } /* for i*/
  for(i=1; i<=NumOfDynVar(); i++) {
    GetDynVarNo(i, name);
    fprintf(pcode, "%s:=xx_[%d];\n", name, i+NoOfSpec()-NumOfConstraint());
  } /* for i */
  for(i=1; i<=NoOfSpec(); i++) {
    GetSpecNo(i, name, &charge);
    RenameSpec(rename, name, charge);
    constraint=IsSpecInConstraint(name, charge);
    if (constraint>0) {
      fprintf(pcode, "%s:=", rename);
      TreePrint(con[constraint-1], 2, pcode);
      fprintf(pcode, ";\n");
    }; /* if */
  }; /* for i */

  for(i=1; i<=NumbOfDynVars; i++) {
    fprintf(pcode, "vv_[%d]:=", i);
    TreePrint(v[i-1], 2, pcode);
    fprintf(pcode, ";\n");
  } /* for i */
  fprintf(pcode, "IF bj_ THEN BEGIN\n");

  for(i=1; i<=NumbOfDynVars; i++) {
    for(j=1; j<=NumbOfDynVars; j++) {
      temp=TreeEval(jacobi[i-1][j-1]);
      if (TreeGetError()==NoEval) {
	fprintf(pcode, "jj_[%d, %d]:=", i, j);
	TreePrint(jacobi[i-1][j-1], 2, pcode);
	fprintf(pcode, ";\n");
      } /* if */
    } /* for j */
  } /* for i */
  fprintf(pcode, "END;\nEND; (* derivs *)\n\n");

  /* End of printing derivs */


  /* Printing derivsinit */

  fprintf(pcode, "PROCEDURE derivsinit;\n");
  fprintf(pcode, "BEGIN\n");

  dyn=1;
  for(i=1; i<=NoOfSpec(); i++) {
    GetSpecNo(i, name, &charge);
    if (IsSpecInConstraint(name, charge)==0) {
      RenameSpec(rename, name, charge);
      fprintf(pcode, "species[%d]:='%s';\n", dyn, rename);
      dyn++;
    } /* if */
  } /* for i */

  /* maaske skal foelgende benyttes istedet for det naeste */
 /* for(i=1; i<=NumOfDynVar(); i++) {
    GetDynVarNo(i, name);
    fprintf(pcode, "species[%d]:='%s';\n",
	    i+NoOfSpec()-NumOfConstraint()-IsNotAutoSystem()-1, name);
  }*/ /* for i */

  /* maaske foelgende godt nok */
  for(i=1; i<=NumOfDynVar(); i++) {
    GetDynVarNo(i, name);
    fprintf(pcode, "species[%d]:='%s';\n", i+NoOfSpec(), name);
  }; /* for i */

  GetAndPrintConst("epsr", "epsr", 1, 1e-3, pcode, 2);
  GetAndPrintConst("epsa", "epsa", 1, 1e-20, pcode, 2);
  GetStrConst("datafile", name);
  if (GetError()==NotFound)
    fprintf(pcode, "name_datafile:='kinwrkdat';\n");
  else
    fprintf(pcode, "name_datafile:='%s';\n", name);
  dyn=1;
  for(i=1; i<=NoOfSpec(); i++) {
    GetSpecNo(i, name, &charge);
    if (IsSpecInConstraint(name, charge)==0) {
      temp=GetBeginConc(name, charge);
      fprintf(pcode, "xx[%d]:=%e;\n", dyn, temp);
      dyn++;
    } /* if */
  } /* for i */
  for(i=1; i<=NumOfDynVar(); i++) {
    GetDynVarNo(i, name);
    fprintf(pcode, "xx[%d]:=%e;\n", i+NoOfSpec()-NumOfConstraint(),
	    GetInitValue(name));
  } /* for i */

  /* constant elements of the Jacobian matrix */
  for(i=1; i<=NumbOfDynVars; i++) {
    for(j=1; j<=NumbOfDynVars; j++) {
      temp=TreeEval(jacobi[i-1][j-1]);
      if (TreeGetError()==TreeNoError) {
	fprintf(pcode, "jacobi[%d, %d]:= %e;\n", i, j, temp);
      } /* if */
    } /* for j */
  } /* for i */

  fprintf(pcode, "END; (* derivs *)\n\n");

  /* End of printing derivsinit */

  /* NumbOfParams: 1 - sp_dalfa and derpinit are printed, 2 - hf_dalfa and hopfinit are printed */

  switch (NumbOfParams) {
  case 1:
    /* Printing sp_dalfa */
    fprintf(pcode, "PROCEDURE sp_dalfa(bj_: BOOLEAN;xx_: glnarray; VAR gg_: glnpbynp);\n");
    fprintf(pcode, "BEGIN\n");
    dyn=1;
    for(i=1; i<=NoOfSpec(); i++) {
      GetSpecNo(i, name, &charge);
      if ((IsSpecInConstraint(name, charge)==0) &&
	  (IsSpecParam(name, charge)!=1)) {
	RenameSpec(rename, name, charge);
	fprintf(pcode, "%s:=xx_[%d];\n", rename, dyn);
	dyn++;
      } /* if */
    } /* for i */
    for(i=1; i<=NumOfDynVar(); i++) {
      GetDynVarNo(i, name);
      fprintf(pcode, "%s:=xx_[%d];\n", name, i+NoOfSpec()-NumOfConstraint());
    } /* for i */

    GetParamNo(1, name, &charge, &form);
    if (form==2)
      RenameSpec(rename, name, charge);
    else
      strcpy(rename, name);
    fprintf(pcode, "%s:=xx_[n1];\n", rename);
    fprintf(pcode, "IF bj_ THEN\nBEGIN\n");
    GetParamNo(1, name, &charge, &form);
    if (form==1)
      strcpy(rename, name);
    else
      RenameSpec(rename, name, charge);
    for(j=0; j<NumbOfDynVars; j++) {
      tmp=TreeCreate();
      TreeDerive(tmp, v[j], rename);
      fprintf(pcode, "gg_[%d, n1]:=", j+1);
      TreePrint(tmp, 2, pcode);
      fprintf(pcode, ";\n");
      TreeKill(tmp);
    } /* for j */
    fprintf(pcode, "END;\nEND; (* sp_dalfa *)\n\n");
    /* End of printing sp_dalfa */

    /* Printing derpinit */

    fprintf(pcode, "PROCEDURE derpinit;\nBEGIN\n");

    GetAndPrintConst("cfout", "cfout", 0, 10.0, pcode, 2);
    GetAndPrintConst("maxnoofp", "maxnoofp", 0, 10.0, pcode, 2);
    GetAndPrintConst("maxithqr2", "maxithqr2", 0, 30.0, pcode, 2);
    GetAndPrintConst("maxititera", "maxititera", 0, 30.0, pcode, 2);
    GetAndPrintConst("maxitcorrec", "maxitcorrec", 0, 30.0, pcode, 2);
    GetAndPrintConst("corrhreg", "corrhreg", 0, 10.0, pcode, 2);
    GetAndPrintConst("hh", "hh", 1, 1.0, pcode, 2);
    GetAndPrintConst("epsfx", "epsfx", 1, 1.0e-15, pcode, 2);
    GetAndPrintConst("epsblk", "epsblk", 1, 1.0e-15, pcode, 2);
    GetAndPrintConst("epshqr2", "epshqr2", 1, 1.0E-20, pcode, 2);
    GetAndPrintConst("epsmach", "epsmach", 1, 1e-15, pcode, 2);
    GetAndPrintConst("epsdigit", "epsdigit", 1, 1e-14, pcode, 2);
    GetAndPrintConst("prnscr", "print_on_screen", 0, 1.0, pcode, 2);
    GetAndPrintConst("prnfile", "print_on_file", 0, 1.0, pcode, 2);
    GetAndPrintConst("failprn", "failure_print", 0, 0.0, pcode, 2);
    GetAndPrintConst("paramregu", "para_regu", 0, 1.0, pcode, 2);
    GetAndPrintConst("vecprn", "eigvec_prn", 0, 0.0, pcode, 2);

    GetStrConst("datafile", name);
    if (GetError()==NotFound)
      fprintf(pcode, "name_datafile:= 'kinwrkdat';\n");
    else
      fprintf(pcode, "name_datafile:= '%s';\n", name);
    GetStrConst("datafile", name);
    if (GetError()==NotFound)
      fprintf(pcode, "name_textfile:= 'kinwrkdat.t';\n");
    else
      fprintf(pcode, "name_textfile:= '%s.t';\n", name);

    dyn=0;
    for(i=1; i<=NoOfSpec(); i++) {
      GetSpecNo(i, name, &charge);
      if (IsSpecInConstraint(name, charge)==0) {
	dyn++;
	fprintf(pcode, "xx[%d]:=%e;\n", dyn, GetBeginConc(name, charge));

	temp=GetSpecConst(name, charge, "Bfstep");
	BfErrorCode= GetError();
	if (BfErrorCode==NotFound) {
	  fprintf(pcode, "initndir[%d]:=-1;\n", dyn);
	  fprintf(pcode, "inithmax[%d]:=%e;\n", dyn,1.0E2);
        } else {
	  if (BfErrorCode==NoError) {
	     if (temp>0.0)
	       fprintf(pcode, "initndir[%d]:=1;\n", dyn);
	     else
	       fprintf(pcode, "initndir[%d]:=-1;\n", dyn);
	     fprintf(pcode, "inithmax[%d]:=%e;\n", dyn, fabs(temp));
	  } /* if */
	} /* else */

	temp=GetSpecConst(name, charge, "Bfxmax");
	BfErrorCode= GetError();
	if (BfErrorCode==NotFound) {
	  fprintf(pcode, "initxupp[%d]:=%e;\n", dyn, 1.0E3);
        } else {
	  if (BfErrorCode==NoError) {
       	     fprintf(pcode, "initxupp[%d]:=%e;\n", dyn, temp);
          } /* if */
        } /* else */

	temp=GetSpecConst(name, charge, "Bfxmin");
	BfErrorCode= GetError();
	if (BfErrorCode==NotFound) {
	  fprintf(pcode, "initxlow[%d]:=%e;\n", dyn, 0.0);
        } else {
	  if (BfErrorCode==NoError) {
	     fprintf(pcode, "initxlow[%d]:=%e;\n", dyn, temp);
          } /* if */
        } /* else */

	temp=GetSpecConst(name, charge, "Bfpref");
	BfErrorCode= GetError();
	if (BfErrorCode==NotFound) {
	  fprintf(pcode, "initpref[%d]:=%e;\n", dyn, 0.1);
        } else {
	  if (BfErrorCode==NoError) {
	     fprintf(pcode, "initpref[%d]:=%e;\n", dyn, temp);
          } /* if */
        } /* else */


      } /* if */
    } /* for i */

    for(i=1; i<=NumOfParameter(); i++) {
      GetParamNo(i, name, &charge, &form);
      if (form==1) {
	GetInitParam(name, &temp);
	fprintf(pcode, "xx[n%d]:=%e;\n", i, temp);
	temp=GetDirectForParam(name);
	if (temp>0.0)
	  fprintf(pcode, "initndir[n%d]:=1;\n", i);
	else
	  fprintf(pcode, "initndir[n%d]:=-1;\n", i);
	GetDeltaParam(name, &temp);
	fprintf(pcode, "inithmax[n%d]:=%e;\n", i, fabs(temp));
	GetLowHighPrefParam(name, &temp, &temp1, &temp2);
	fprintf(pcode, "initxlow[n%d]:=%e;\n", i, temp);
	fprintf(pcode, "initxupp[n%d]:=%e;\n", i, temp1);
	if ((temp2>1.0) && (temp2<0.0))
	  fprintf(pcode, "KnCont: Pref. for %s must between 0 and 1.\n", name);
	else
	  fprintf(pcode, "initpref[n%d]:=%e;\n", i, temp2);
      } /* if */
      else {
	GetInitConc(name, charge, &temp);
	fprintf(pcode, "xx[n%d]:=%e;\n", i, temp);
	temp=GetDirectForConc(name, charge);
	if (temp>0.0)
	  fprintf(pcode, "initndir[n%d]:=1;\n", i);
	else
	  fprintf(pcode, "initndir[n%d]:=-1;\n", i);
	GetDeltaConc(name, charge, &temp);
	fprintf(pcode, "inithmax[n%d]:=%e;\n", i, fabs(temp));
	GetLowHighPrefConc(name, charge, &temp, &temp1, &temp2);
	fprintf(pcode, "initxlow[n%d]:=%e;\n", i, temp);
	fprintf(pcode, "initxupp[n%d]:=%e;\n", i, temp1);
	if ((temp2>1.0) && (temp2<0.0))
	  fprintf(pcode, "KnCont: Pref. for %s(%d) must between 0 and 1.\n",
		  name, (int)charge);
	else
	  fprintf(pcode, "initpref[n%d]:=%e;\n", i, temp2);
      } /* else */
    } /* for i */
    fprintf(pcode, "END; (* derpinit *)\n\n");
    /* End of printing derpinit */

    fprintf(pcode, "PROCEDURE hf_dalfa(bj_,tp_: BOOLEAN;xx_: glnarray; VAR gg_: glnpbynp);\n");
    fprintf(pcode, "BEGIN\n");
    fprintf(pcode, "END; (* hf_dalfa *)\n\n");
    /* End of printing hf_dalfa */

    /*Printing hopfinit */
    fprintf(pcode, "PROCEDURE hopfinit;\nBEGIN\n");
    fprintf(pcode, "END; (* hopfinit *)\n\n");
    /* End of printing hopfinit */

    break;
  case 2:
    /*Printing hf_dalfa */

    fprintf(pcode, "PROCEDURE sp_dalfa(bj_: BOOLEAN;xx_: glnarray; VAR gg_: glnpbynp);\n");
    fprintf(pcode, "BEGIN\n");
    fprintf(pcode, "END; (* sp_dalfa *)\n\n");
    /* End of printing sp_dalfa */

    /* Printing derpinit */
    fprintf(pcode, "PROCEDURE derpinit;\nBEGIN\n");
    fprintf(pcode, "END; (* derpinit *)\n\n");
    /* End of printing derpinit */

    fprintf(pcode, "PROCEDURE hf_dalfa(bj_,tp_: BOOLEAN;xx_: glnarray; VAR gg_: glnpbynp);\n");
    fprintf(pcode, "BEGIN\n");
    dyn=1;
    for(i=1; i<=NoOfSpec(); i++) {
      GetSpecNo(i, name, &charge);
      if ((IsSpecInConstraint(name, charge)==0) &&
	  (IsSpecParam(name, charge)!=1)) {
	RenameSpec(rename, name, charge);
	fprintf(pcode, "%s:=xx_[%d];\n", rename, dyn);
	dyn++;
      } /* if */
    } /* for i */
    for(i=1; i<=NumOfDynVar(); i++) {
      GetDynVarNo(i, name);
      fprintf(pcode, "%s:=xx_[%d];\n", name, i+NoOfSpec()-NumOfConstraint());
    } /* for i */

    fprintf(pcode, "IF tp_ THEN BEGIN\n");
    for(i=1; i<=NumbOfParams; i++) {
      GetParamNo(i, name, &charge, &form);
      if (form==2)
	RenameSpec(rename, name, charge);
      else
	strcpy(rename, name);
      fprintf(pcode, "%s:=xx_[n%d];\n", rename, i);
    } /* for i */
    fprintf(pcode, "IF bj_ THEN BEGIN\n");
    GetParamNo(1, name, &charge, &form);
    if (form==2)
      RenameSpec(rename, name, charge);
    else
      strcpy(rename, name);
    for(j=0; j<NumbOfDynVars; j++) {
      tmp=TreeCreate();
      TreeDerive(tmp, v[j], rename);
      fprintf(pcode, "gg_[%d, n1]:=", j+1);
      TreePrint(tmp, 2, pcode);
      fprintf(pcode, ";\n");
      TreeKill(tmp);
    } /* for j */
    fprintf(pcode, "END;\nEND\nELSE BEGIN\n");
    for(i=NumbOfParams; i>=1; i--) {
      GetParamNo(i, name, &charge, &form);
      if (form==2)
	RenameSpec(rename, name, charge);
      else
	strcpy(rename, name);
      fprintf(pcode, "%s:=xx_[n%d];\n", rename, NumbOfParams+1-i);
    } /* for i */
    fprintf(pcode, "IF bj_ THEN BEGIN\n");
    GetParamNo(2, name, &charge, &form);
    if (form==2)
      RenameSpec(rename, name, charge);
    else
      strcpy(rename, name);
    for(j=0; j<NumbOfDynVars; j++) {
      tmp=TreeCreate();
      TreeDerive(tmp, v[j], rename);
      fprintf(pcode, "gg_[%d, n1]:=", j+1);
      TreePrint(tmp, 2, pcode);
      fprintf(pcode, ";\n");
      TreeKill(tmp);
    } /* for i */
    fprintf(pcode, "END;\nEND;\nEND; (* hf_dalfa *)\n\n");

    /* End of printing hf_dalfa */

    /*Printing hopfinit */

    fprintf(pcode, "PROCEDURE hopfinit;\nBEGIN\n");
    fprintf(pcode, "need_dd_jac:=%s;\n",
	    (need_dd_jac==1)?"TRUE":"FALSE");

    GetAndPrintConst("maxoutn1", "maxoutn1", 0, 30.0, pcode, 2);
    GetAndPrintConst("maxoutn2", "maxoutn2", 0, 2.0, pcode, 2);
    GetAndPrintConst("cfout", "cfout", 0, 10.0, pcode, 2);
    GetAndPrintConst("ref", "ref", 0, 1.0, pcode, 2);
    GetAndPrintConst("hopfbftp", "hopfbftp", 0, -1.0, pcode, 2);
    GetAndPrintConst("re1", "re1", 1, 0.5, pcode, 2);
    GetAndPrintConst("im1", "im1", 1, 0.5, pcode, 2);
    GetAndPrintConst("maxnoofp", "maxnoofp", 0, 10.0, pcode, 2);
    GetAndPrintConst("maxithqr2", "maxithqr2", 0, 30.0, pcode, 2);
    GetAndPrintConst("maxitbisec", "maxitbisec", 0, 30.0, pcode, 2);
    GetAndPrintConst("maxitintpol", "maxitintpol", 0, 30.0, pcode, 2);
    GetAndPrintConst("maxititera", "maxititera", 0, 30.0, pcode, 2);
    GetAndPrintConst("maxitcorrec", "maxitcorrec", 0, 30.0, pcode, 2);
    GetAndPrintConst("corrhreg", "corrhreg", 0, 10.0, pcode, 2);
    GetAndPrintConst("hh", "hh", 1, 1.0, pcode, 2);
    GetAndPrintConst("epsfx", "epsfx", 1, 1.0e-15, pcode, 2);
    GetAndPrintConst("epsblk", "epsblk", 1, 1.0e-15, pcode, 2);
    GetAndPrintConst("epshqr2", "epshqr2", 1, 1.0E-20, pcode, 2);
    GetAndPrintConst("epsbisecr", "epsbisecr", 1, 1.0e-10, pcode, 2);
    GetAndPrintConst("epsmach", "epsmach", 1, 1e-15, pcode, 2);
    GetAndPrintConst("epsdigit", "epsdigit", 1, 1e-14, pcode, 2);
    GetAndPrintConst("prnscr", "print_on_screen", 0, 1.0, pcode, 2);
    GetAndPrintConst("prnfile", "print_on_file", 0, 1.0, pcode, 2);
    GetAndPrintConst("failprn", "failure_print", 0, 0.0, pcode, 2);
    GetAndPrintConst("paramregu", "para_regU", 0, 1.0, pcode, 2);
    GetAndPrintConst("qccalc", "qc_calc", 0, 0.0, pcode, 2);
    GetAndPrintConst("hfcalc", "hf_calc", 0, 1.0, pcode, 2);
    GetAndPrintConst("dfcalc", "df_calc", 0, 0.0, pcode, 2);
    GetAndPrintConst("hasscalc", "hass_calc", 0, 1.0, pcode, 2);
    GetAndPrintConst("kuracalc", "kura_calc", 0, 0.0, pcode, 2);
    GetAndPrintConst("vecprn", "eigvec_prn", 0, 0.0, pcode, 2);

    GetStrConst("datafile", name);
    if (GetError()==NotFound)
      fprintf(pcode, "name_datafile:= 'kinwrkdat';\n");
    else
      fprintf(pcode, "name_datafile:= '%s';\n", name);
    GetStrConst("datafile", name);
    if (GetError()==NotFound)
      fprintf(pcode, "name_textfile:= 'kinwrkdat.t';\n");
    else
      fprintf(pcode, "name_textfile:= '%s.t';\n", name);

    dyn=0;
    for(i=1; i<=NoOfSpec(); i++) {
      GetSpecNo(i, name, &charge);
      if (IsSpecInConstraint(name, charge)==0) {
	dyn++;
	fprintf(pcode, "xx[%d]:=%e;\n", dyn, GetBeginConc(name, charge));

	temp=GetSpecConst(name, charge, "Bfstep");
	BfErrorCode= GetError();
	if (BfErrorCode==NotFound) {
	  fprintf(pcode, "initndir[%d]:=-1;\n", dyn);
	  fprintf(pcode, "inithmax[%d]:=%e;\n", dyn,1.0E2);
        } else {
	  if (BfErrorCode==NoError) {
	     if (temp>0.0)
	       fprintf(pcode, "initndir[%d]:=1;\n", dyn);
	     else
	       fprintf(pcode, "initndir[%d]:=-1;\n", dyn);
	     fprintf(pcode, "inithmax[%d]:=%e;\n", dyn, fabs(temp));
	  } /* if */
	} /* else */

	temp=GetSpecConst(name, charge, "Bfxmax");
	BfErrorCode= GetError();
	if (BfErrorCode==NotFound) {
	  fprintf(pcode, "initxupp[%d]:=%e;\n", dyn, 1.0E3);
        } else {
	  if (BfErrorCode==NoError) {
       	     fprintf(pcode, "initxupp[%d]:=%e;\n", dyn, temp);
          } /* if */
        } /* else */

	temp=GetSpecConst(name, charge, "Bfxmin");
	BfErrorCode= GetError();
	if (BfErrorCode==NotFound) {
	  fprintf(pcode, "initxlow[%d]:=%e;\n", dyn, 0.0);
        } else {
	  if (BfErrorCode==NoError) {
	     fprintf(pcode, "initxlow[%d]:=%e;\n", dyn, temp);
          } /* if */
        } /* else */

	temp=GetSpecConst(name, charge, "Bfpref");
	BfErrorCode= GetError();
	if (BfErrorCode==NotFound) {
	  fprintf(pcode, "initpref[%d]:=%e;\n", dyn, 0.1);
        } else {
	  if (BfErrorCode==NoError) {
	     fprintf(pcode, "initpref[%d]:=%e;\n", dyn, temp);
          } /* if */
        } /* else */

      } /* if */
    } /* for i */
    for(i=1; i<=NumOfParameter(); i++) {
      GetParamNo(i, name, &charge, &form);
      if (form==1) {
	GetInitParam(name, &temp);
	fprintf(pcode, "xx[n%d]:=%e;\n", i, temp);
	temp=GetDirectForParam(name);
	if (temp>0.0)
	  fprintf(pcode, "initndir[n%d]:=1;\n", i);
	else
	  fprintf(pcode, "initndir[n%d]:=-1;\n", i);
	GetDeltaParam(name, &temp);
	fprintf(pcode, "inithmax[n%d]:=%e;\n", i, fabs(temp));
	GetLowHighPrefParam(name, &temp, &temp1, &temp2);
	fprintf(pcode, "initxlow[n%d]:=%e;\n", i, temp);
	fprintf(pcode, "initxupp[n%d]:=%e;\n", i, temp1);
	if ((temp2>1.0) && (temp2<0.0))
	  fprintf(pcode, "KnCont: Pref. for %s must between 0 and 1.\n", name);
	else
	  fprintf(pcode, "initpref[n%d]:=%e;\n", i, temp2);
      } /* if */
      else {
	GetInitConc(name, charge, &temp);
	fprintf(pcode, "xx[n%d]:=%e;\n", i, temp);
	temp=GetDirectForConc(name, charge);
	if (temp>0.0)
	  fprintf(pcode, "initndir[n%d]:=1;\n", i);
	else
	  fprintf(pcode, "initndir[n%d]:=-1;\n", i);
	GetDeltaConc(name, charge, &temp);
	fprintf(pcode, "inithmax[n%d]:=%e;\n", i, fabs(temp));
	GetLowHighPrefConc(name, charge, &temp, &temp1, &temp2);
	fprintf(pcode, "initxlow[n%d]:=%e;\n", i, temp);
	fprintf(pcode, "initxupp[n%d]:=%e;\n", i, temp1);
	if ((temp2>1.0) && (temp2<0.0))
	  fprintf(pcode, "KnCont: Pref. for %s(%d) must between 0 and 1.\n",
		  name, (int)charge);
	else
	  fprintf(pcode, "initpref[n%d]:=%e;\n", i, temp2);
      } /* else */
    } /* for i */
    fprintf(pcode, "END; (* hopfinit *)\n\n");

    /* End of printing hopfinit */

    break;
  } /* switch NumbOfParams */

  /* Printing detnumparam */

  fprintf(pcode, "PROCEDURE detnumparam;\n");
  fprintf(pcode, "BEGIN\n");
  fprintf(pcode, "  numparam:=%d;\n", NumbOfParams);
  fprintf(pcode, "END; (* detnumparam *)\n\n");
  /* End of printing detnumparam */


  StringFree(name);
  StringFree(rename);
} /* KnCont */
