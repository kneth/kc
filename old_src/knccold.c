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

  fprintf(ccode, "(******************************************************\n");
  fprintf(ccode, "  WARNING: This file was generated by kc v%s\n",
	  VERSION);
  fprintf(ccode, "  CopyWrong 1994 by Kenneth Geisshirt.\n");
  fprintf(ccode, "  %s", ctime(&timer));
  fprintf(ccode, "*******************************************************)\n");

  /* printing Hessian */
  fprintf(ccode, "void djacobian(double *xx_, double ***dS_)\n");
  fprintf(ccode, "{\n");
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
	fprintf(ccode, "%s=x(%d);\n", rename, dyn);
	dyn++;
      } /* if */
    } /* for i */
    for(i=1; i<=NumOfDynVar(); i++) {
      GetDynVarNo(i, name);
      fprintf(ccode, "%s=x(%d);\n", rename,
	      i+NoOfSpec()/NumOfConstraint());
    } /* for i */
    need_dd_jac=0;

/*
    for(i=0; i<NumbOfDynVars; i++) {
      fprintf(ccode, "dS_[%d,%d,%d]=", i+1, i+1, i+1);
      temp=TreeEval(hess[i][i][i]);
      if (TreeGetError()==NoEval) {
	need_dd_jac=1;
	TreePrint(hess[i][i][i], 3, ccode);
      } else
	fprintf(ccode, "%e", temp);
      fprintf(ccode, ";\n");
    }*/ /* for i */

    for(i=0; i<NumbOfDynVars; i++)
      for(j=0; j<NumbOfDynVars; j++)
      /*
	for(l=0; l<=j; l++) {
	*/
	for(l=0; l<NumbOfDynVars; l++) {
	  fprintf(ccode, "dS_[%d,%d,%d]=", i+1, j+1, l+1);
	  temp=TreeEval(hess[i][j][l]);
	  if (TreeGetError()==NoEval) {
	    need_dd_jac=1;
	    TreePrint(hess[i][j][l], 3, ccode);
	  } else
	    fprintf(ccode, "%e", temp);
	  fprintf(ccode, ";\n");
	  /*
	  if (l<j) {
	  fprintf(ccode, "dS_[%d,%d,%d]=dS_[%d,%d,%d];\n", i+1, l+1,
		  j+1, i+1, j+1, l+1);
          }*/ /* if */
	} /* for l */
    break;
  } /* switch NumbOfParams */
  fprintf(ccode, "} /* djacobian */\n\n");

  /* End of printing the Hessian */

  /* printing Keldian */
  fprintf(ccode, "void ddjacobian(double *xx_,double **** ddS_)\n");
  fprintf(ccode, "{\n");
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
      fprintf(ccode, "%s=x(%d);\n", rename, dyn);
      dyn++;
    } /* if */
  } /* for i */
  for(i=1; i<=NumOfDynVar(); i++) {
    GetDynVarNo(i, name);
    fprintf(ccode, "%s=x(%d);\n", name, i+NoOfSpec()-NumOfConstraint());
  } /* for i */
    /*
    for(i=0; i<NumbOfDynVars; i++) {
      fprintf(ccode, "dds_[%d,%d,%d,%d]=", i+1, i+1, i+1, i+1);
      TreePrint(keld[i][i][i][i], 3, ccode);
      fprintf(ccode, ";\n");
    }*/ /* i */


    for(i=0; i<NumbOfDynVars; i++)
      for(j=0; j<NumbOfDynVars; j++)
	for(l=0; l<NumbOfDynVars; l++)
	/*
	  for(k=0; k<=l; k++) {
	  */
	  for(k=0; k<NumbOfDynVars; k++) {
	    fprintf(ccode, "dds_[%d,%d,%d,%d]=", i+1, j+1, l+1, k+1);
	    TreePrint(keld[i][j][l][k], 3, ccode);
	    fprintf(ccode, ";\n");
	    /*
	    fprintf(ccode, "dds_[%d,%d,%d,%d]=ddS_[%d,%d,%d,%d];\n",
		    i+1, j+1, k+1, l+1, i+1, j+1, l+1, k+1);
	    fprintf(ccode, "dds_[%d,%d,%d,%d]=ddS_[%d,%d,%d,%d];\n",
		    i+1, l+1, j+1, k+1, i+1, j+1, l+1, k+1);
	    fprintf(ccode, "dds_[%d,%d,%d,%d]=ddS_[%d,%d,%d,%d];\n",
		    i+1, l+1, k+1, j+1, i+1, j+1, l+1, k+1);
	    fprintf(ccode, "dds_[%d,%d,%d,%d]=ddS_[%d,%d,%d,%d];\n",
		    i+1, k+1, l+1, j+1, i+1, j+1, l+1, k+1);
	    fprintf(ccode, "dds_[%d,%d,%d,%d]=ddS_[%d,%d,%d,%d];\n",
		    i+1, k+1, j+1, l+1, i+1, j+1, l+1, k+1);
	  */
	  } /* for k */
  } /* need_dd_jac==1 */
  fprintf(ccode, "} /* ddjacobian */\n\n");

  /* End of printing the Keldian */


  /* printing derivs */

  fprintf(ccode, "void derivs(int bj_, double xx_, double t_,\n");
  fprintf(ccode, "                 double *vv_, double **jj_)\n");
  fprintf(ccode, "{\n");
  dyn=1;
  for(i=1; i<=NoOfSpec(); i++) {
    GetSpecNo(i, name, &charge);
    RenameSpec(rename, name, charge);
    constraint=IsSpecInConstraint(name, charge);
    if (constraint==0) {
      fprintf(ccode, "%s=x(%d);\n", rename, dyn);
      dyn++;
    } /* if */
  } /* for i*/
  for(i=1; i<=NumOfDynVar(); i++) {
    GetDynVarNo(i, name);
    fprintf(ccode, "%s=x(%d);\n", name, i+NoOfSpec()-NumOfConstraint());
  } /* for i */
  for(i=1; i<=NoOfSpec(); i++) {
    GetSpecNo(i, name, &charge);
    RenameSpec(rename, name, charge);
    constraint=IsSpecInConstraint(name, charge);
    if (constraint>0) {
      fprintf(ccode, "%s=", rename);
      TreePrint(con[constraint-1], 3, ccode);
      fprintf(ccode, ";\n");
    }; /* if */
  }; /* for i */

  for(i=1; i<=NumbOfDynVars; i++) {
    fprintf(ccode, "f(%d) = ", i);
    TreePrint(v[i-1], 3, ccode);
    fprintf(ccode, ";\n");
  } /* for i */
  fprintf(ccode, "if (bj_==1) {\n");

  for(i=1; i<=NumbOfDynVars; i++) {
    for(j=1; j<=NumbOfDynVars; j++) {
      temp=TreeEval(jacobi[i-1][j-1]);
  /*    if (TreeGetError()==NoEval) {
*/	fprintf(ccode, "g(%d,%d) = ", i, j);
	TreePrint(jacobi[i-1][j-1], 3, ccode);
	fprintf(ccode, ";\n");
   /*   } /* if */
   
    } /* for j */
  } /* for i */
  fprintf(ccode, "}\n} /* derivs */\n\n");

  /* End of printing derivs */


  /* Printing derivsinit */

  fprintf(ccode, "void derivsinit(void)\n");
  fprintf(ccode, "{\n");

  dyn=1;
  for(i=1; i<=NoOfSpec(); i++) {
    GetSpecNo(i, name, &charge);
    if (IsSpecInConstraint(name, charge)==0) {
      RenameSpec(rename, name, charge);
      fprintf(ccode, "species[%d]='%s';\n", dyn, rename);
      dyn++;
    } /* if */
  } /* for i */

  /* maaske skal foelgende benyttes istedet for det naeste */
 /* for(i=1; i<=NumOfDynVar(); i++) {
    GetDynVarNo(i, name);
    fprintf(ccode, "species[%d]='%s';\n",
	    i+NoOfSpec()-NumOfConstraint()-IsNotAutoSystem()-1, name);
  }*/ /* for i */

  /* maaske foelgende godt nok */
  for(i=1; i<=NumOfDynVar(); i++) {
    GetDynVarNo(i, name);
    fprintf(ccode, "species[%d]='%s';\n", i+NoOfSpec(), name);
  }; /* for i */

  GetAndPrintConst("epsr", "epsr", 1, 1e-3, ccode, 3);
  GetAndPrintConst("epsa", "epsa", 1, 1e-20, ccode, 3);
  GetStrConst("datafile", name);
  if (GetError()==NotFound)
    fprintf(ccode, "name_datafile='kinwrkdat';\n");
  else
    fprintf(ccode, "name_datafile='%s';\n", name);
  dyn=1;
  for(i=1; i<=NoOfSpec(); i++) {
    GetSpecNo(i, name, &charge);
    if (IsSpecInConstraint(name, charge)==0) {
      temp=GetBeginConc(name, charge);
      fprintf(ccode, "xx[%d]=%e;\n", dyn, temp);
      dyn++;
    } /* if */
  } /* for i */
  for(i=1; i<=NumOfDynVar(); i++) {
    GetDynVarNo(i, name);
    fprintf(ccode, "xx[%d]=%e;\n", i+NoOfSpec()-NumOfConstraint(),
	    GetInitValue(name));
  } /* for i */

  /* constant elements of the Jacobian matrix */
  for(i=1; i<=NumbOfDynVars; i++) {
    for(j=1; j<=NumbOfDynVars; j++) {
      temp=TreeEval(jacobi[i-1][j-1]);
      if (TreeGetError()==TreeNoError) {
	fprintf(ccode, "jacobi[%d, %d]= %e;\n", i, j, temp);
      } /* if */
    } /* for j */
  } /* for i */

  fprintf(ccode, "} /* derivs */\n\n");

  /* End of printing derivsinit */

  /* NumbOfParams: 1 - sp_dalfa and derpinit are printed, 2 - hf_dalfa and hopfinit are printed */

  switch (NumbOfParams) {
  case 1:
    /* Printing sp_dalfa */
    fprintf(ccode, "void sp_dalfa(int bj_,double xx_, double **gg_)\n");
    fprintf(ccode, "{\n");
    dyn=1;
    for(i=1; i<=NoOfSpec(); i++) {
      GetSpecNo(i, name, &charge);
      if ((IsSpecInConstraint(name, charge)==0) &&
	  (IsSpecParam(name, charge)!=1)) {
	RenameSpec(rename, name, charge);
	fprintf(ccode, "%s=xx_[%d];\n", rename, dyn);
	dyn++;
      } /* if */
    } /* for i */
    for(i=1; i<=NumOfDynVar(); i++) {
      GetDynVarNo(i, name);
      fprintf(ccode, "%s=xx_[%d];\n", name, i+NoOfSpec()-NumOfConstraint());
    } /* for i */

    GetParamNo(1, name, &charge, &form);
    if (form==2)
      RenameSpec(rename, name, charge);
    else
      strcpy(rename, name);
    fprintf(ccode, "%s=xx_[n1];\n", rename);
    fprintf(ccode, "if (bj_==1) \n{\n");
    GetParamNo(1, name, &charge, &form);
    if (form==1)
      strcpy(rename, name);
    else
      RenameSpec(rename, name, charge);
    for(j=0; j<NumbOfDynVars; j++) {
      tmp=TreeCreate();
      TreeDerive(tmp, v[j], rename);
      fprintf(ccode, "gg_[%d, n1]=", j+1);
      TreePrint(tmp, 3, ccode);
      fprintf(ccode, ";\n");
      TreeKill(tmp);
    } /* for j */
    fprintf(ccode, "}\n} /* sp_dalfa */\n\n");
    /* End of printing sp_dalfa */

    /* Printing derpinit */

    fprintf(ccode, "void derpinit(void)\n{\n");

    GetAndPrintConst("cfout", "cfout", 0, 10.0, ccode, 3);
    GetAndPrintConst("maxnoofp", "maxnoofp", 0, 10.0, ccode, 3);
    GetAndPrintConst("maxithqr2", "maxithqr2", 0, 30.0, ccode, 3);
    GetAndPrintConst("maxititera", "maxititera", 0, 30.0, ccode, 3);
    GetAndPrintConst("maxitcorrec", "maxitcorrec", 0, 30.0, ccode, 3);
    GetAndPrintConst("corrhreg", "corrhreg", 0, 10.0, ccode, 3);
    GetAndPrintConst("hh", "hh", 1, 1.0, ccode, 3);
    GetAndPrintConst("epsfx", "epsfx", 1, 1.0e-15, ccode, 3);
    GetAndPrintConst("epsblk", "epsblk", 1, 1.0e-15, ccode, 3);
    GetAndPrintConst("epshqr2", "epshqr2", 1, 1.0E-20, ccode, 3);
    GetAndPrintConst("epsmach", "epsmach", 1, 1e-15, ccode, 3);
    GetAndPrintConst("epsdigit", "epsdigit", 1, 1e-14, ccode, 3);
    GetAndPrintConst("prnscr", "print_on_screen", 0, 1.0, ccode, 3);
    GetAndPrintConst("prnfile", "print_on_file", 0, 1.0, ccode, 3);
    GetAndPrintConst("failprn", "failure_print", 0, 0.0, ccode, 3);
    GetAndPrintConst("paramregu", "para_regu", 0, 1.0, ccode, 3);
    GetAndPrintConst("vecprn", "eigvec_prn", 0, 0.0, ccode, 3);

    GetStrConst("datafile", name);
    if (GetError()==NotFound)
      fprintf(ccode, "name_datafile= 'kinwrkdat';\n");
    else
      fprintf(ccode, "name_datafile= '%s';\n", name);
    GetStrConst("datafile", name);
    if (GetError()==NotFound)
      fprintf(ccode, "name_textfile= 'kinwrkdat.t';\n");
    else
      fprintf(ccode, "name_textfile= '%s.t';\n", name);

    dyn=0;
    for(i=1; i<=NoOfSpec(); i++) {
      GetSpecNo(i, name, &charge);
      if (IsSpecInConstraint(name, charge)==0) {
	dyn++;
	fprintf(ccode, "xx[%d]=%e;\n", dyn, GetBeginConc(name, charge));

	temp=GetSpecConst(name, charge, "Bfstep");
	BfErrorCode= GetError();
	if (BfErrorCode==NotFound) {
	  fprintf(ccode, "initndir[%d]=-1;\n", dyn);
	  fprintf(ccode, "inithmax[%d]=%e;\n", dyn,1.0E2);
        } else {
	  if (BfErrorCode==NoError) {
	     if (temp>0.0)
	       fprintf(ccode, "initndir[%d]=1;\n", dyn);
	     else
	       fprintf(ccode, "initndir[%d]=-1;\n", dyn);
	     fprintf(ccode, "inithmax[%d]=%e;\n", dyn, fabs(temp));
	  } /* if */
	} /* else */

	temp=GetSpecConst(name, charge, "Bfxmax");
	BfErrorCode= GetError();
	if (BfErrorCode==NotFound) {
	  fprintf(ccode, "initxupp[%d]=%e;\n", dyn, 1.0E3);
        } else {
	  if (BfErrorCode==NoError) {
       	     fprintf(ccode, "initxupp[%d]=%e;\n", dyn, temp);
          } /* if */
        } /* else */

	temp=GetSpecConst(name, charge, "Bfxmin");
	BfErrorCode= GetError();
	if (BfErrorCode==NotFound) {
	  fprintf(ccode, "initxlow[%d]=%e;\n", dyn, 0.0);
        } else {
	  if (BfErrorCode==NoError) {
	     fprintf(ccode, "initxlow[%d]=%e;\n", dyn, temp);
          } /* if */
        } /* else */

	temp=GetSpecConst(name, charge, "Bfpref");
	BfErrorCode= GetError();
	if (BfErrorCode==NotFound) {
	  fprintf(ccode, "initpref[%d]=%e;\n", dyn, 0.1);
        } else {
	  if (BfErrorCode==NoError) {
	     fprintf(ccode, "initpref[%d]=%e;\n", dyn, temp);
          } /* if */
        } /* else */


      } /* if */
    } /* for i */

    for(i=1; i<=NumOfParameter(); i++) {
      GetParamNo(i, name, &charge, &form);
      if (form==1) {
	GetInitParam(name, &temp);
	fprintf(ccode, "xx[n%d]=%e;\n", i, temp);
	temp=GetDirectForParam(name);
	if (temp>0.0)
	  fprintf(ccode, "initndir[n%d]=1;\n", i);
	else
	  fprintf(ccode, "initndir[n%d]=-1;\n", i);
	GetDeltaParam(name, &temp);
	fprintf(ccode, "inithmax[n%d]=%e;\n", i, fabs(temp));
	GetLowHighPrefParam(name, &temp, &temp1, &temp2);
	fprintf(ccode, "initxlow[n%d]=%e;\n", i, temp);
	fprintf(ccode, "initxupp[n%d]=%e;\n", i, temp1);
	if ((temp2>1.0) && (temp2<0.0))
	  fprintf(ccode, "KnCont: Pref. for %s must between 0 and 1.\n", name);
	else
	  fprintf(ccode, "initpref[n%d]=%e;\n", i, temp2);
      } /* if */
      else {
	GetInitConc(name, charge, &temp);
	fprintf(ccode, "xx[n%d]=%e;\n", i, temp);
	temp=GetDirectForConc(name, charge);
	if (temp>0.0)
	  fprintf(ccode, "initndir[n%d]=1;\n", i);
	else
	  fprintf(ccode, "initndir[n%d]=-1;\n", i);
	GetDeltaConc(name, charge, &temp);
	fprintf(ccode, "inithmax[n%d]=%e;\n", i, fabs(temp));
	GetLowHighPrefConc(name, charge, &temp, &temp1, &temp2);
	fprintf(ccode, "initxlow[n%d]=%e;\n", i, temp);
	fprintf(ccode, "initxupp[n%d]=%e;\n", i, temp1);
	if ((temp2>1.0) && (temp2<0.0))
	  fprintf(ccode, "KnCont: Pref. for %s(%d) must between 0 and 1.\n",
		  name, (int)charge);
	else
	  fprintf(ccode, "initpref[n%d]=%e;\n", i, temp2);
      } /* else */
    } /* for i */
    fprintf(ccode, "} /* derpinit */n\n");
    /* End of printing derpinit */

    fprintf(ccode, "void hf_dalfa(int bj_,int tp_, double xx_, double **gg_)\n");
    fprintf(ccode, "{\n");
    fprintf(ccode, "} /* hf_dalfa */\n\n");
    /* End of printing hf_dalfa */

    /*Printing hopfinit */
    fprintf(ccode, "void hopfinit(void)\n{\n");
    fprintf(ccode, "} /* hopfinit */\n\n");
    /* End of printing hopfinit */

    break;
  case 2:
    /*Printing hf_dalfa */

    fprintf(ccode, "void sp_dalfa(int bj_,double xx_, double **gg_)\n");
    fprintf(ccode, "{\n");
    fprintf(ccode, "} /* sp_dalfa */\n\n");
    /* End of printing sp_dalfa */

    /* Printing derpinit */
    fprintf(ccode, "void derpinit(void)\n{\n");
    fprintf(ccode, "} /* derpinit */\n\n");
    /* End of printing derpinit */

    fprintf(ccode, "void hf_dalfa(int bj_,int tp_, double xx_, double **gg_: glnpbynp)\n");
    fprintf(ccode, "{\n");
    dyn=1;
    for(i=1; i<=NoOfSpec(); i++) {
      GetSpecNo(i, name, &charge);
      if ((IsSpecInConstraint(name, charge)==0) &&
	  (IsSpecParam(name, charge)!=1)) {
	RenameSpec(rename, name, charge);
	fprintf(ccode, "%s=xx_[%d];\n", rename, dyn);
	dyn++;
      } /* if */
    } /* for i */
    for(i=1; i<=NumOfDynVar(); i++) {
      GetDynVarNo(i, name);
      fprintf(ccode, "%s=xx_[%d];\n", name, i+NoOfSpec()-NumOfConstraint());
    } /* for i */

    fprintf(ccode, "if (tp_==1) {\n");
    for(i=1; i<=NumbOfParams; i++) {
      GetParamNo(i, name, &charge, &form);
      if (form==2)
	RenameSpec(rename, name, charge);
      else
	strcpy(rename, name);
      fprintf(ccode, "%s=xx_[n%d];\n", rename, i);
    } /* for i */
    fprintf(ccode, "if (bj_==1) {\n");
    GetParamNo(1, name, &charge, &form);
    if (form==2)
      RenameSpec(rename, name, charge);
    else
      strcpy(rename, name);
    for(j=0; j<NumbOfDynVars; j++) {
      tmp=TreeCreate();
      TreeDerive(tmp, v[j], rename);
      fprintf(ccode, "gg_[%d, n1]=", j+1);
      TreePrint(tmp, 3, ccode);
      fprintf(ccode, ";\n");
      TreeKill(tmp);
    } /* for j */
    fprintf(ccode, "}\n}\nelse {\n");
    for(i=NumbOfParams; i>=1; i--) {
      GetParamNo(i, name, &charge, &form);
      if (form==2)
	RenameSpec(rename, name, charge);
      else
	strcpy(rename, name);
      fprintf(ccode, "%s=xx_[n%d];\n", rename, NumbOfParams+1-i);
    } /* for i */
    fprintf(ccode, "if (bj_==1) {\n");
    GetParamNo(2, name, &charge, &form);
    if (form==2)
      RenameSpec(rename, name, charge);
    else
      strcpy(rename, name);
    for(j=0; j<NumbOfDynVars; j++) {
      tmp=TreeCreate();
      TreeDerive(tmp, v[j], rename);
      fprintf(ccode, "gg_[%d, n1]=", j+1);
      TreePrint(tmp, 3, ccode);
      fprintf(ccode, ";\n");
      TreeKill(tmp);
    } /* for i */
    fprintf(ccode, "}\n}\n} /* hf_dalfa */\n\n");

    /* End of printing hf_dalfa */

    /*Printing hopfinit */

    fprintf(ccode, "void hopfinit(void)\n{\n");
    fprintf(ccode, "need_dd_jac=%s;\n",
	    (need_dd_jac==1)?"TRUE":"FALSE");

    GetAndPrintConst("maxoutn1", "maxoutn1", 0, 30.0, ccode, 3);
    GetAndPrintConst("maxoutn2", "maxoutn2", 0, 2.0, ccode, 3);
    GetAndPrintConst("cfout", "cfout", 0, 10.0, ccode, 3);
    GetAndPrintConst("ref", "ref", 0, 1.0, ccode, 3);
    GetAndPrintConst("hopfbftp", "hopfbftp", 0, -1.0, ccode, 3);
    GetAndPrintConst("re1", "re1", 1, 0.5, ccode, 3);
    GetAndPrintConst("im1", "im1", 1, 0.5, ccode, 3);
    GetAndPrintConst("maxnoofp", "maxnoofp", 0, 10.0, ccode, 3);
    GetAndPrintConst("maxithqr2", "maxithqr2", 0, 30.0, ccode, 3);
    GetAndPrintConst("maxitbisec", "maxitbisec", 0, 30.0, ccode, 3);
    GetAndPrintConst("maxitintpol", "maxitintpol", 0, 30.0, ccode, 3);
    GetAndPrintConst("maxititera", "maxititera", 0, 30.0, ccode, 3);
    GetAndPrintConst("maxitcorrec", "maxitcorrec", 0, 30.0, ccode, 3);
    GetAndPrintConst("corrhreg", "corrhreg", 0, 10.0, ccode, 3);
    GetAndPrintConst("hh", "hh", 1, 1.0, ccode, 3);
    GetAndPrintConst("epsfx", "epsfx", 1, 1.0e-15, ccode, 3);
    GetAndPrintConst("epsblk", "epsblk", 1, 1.0e-15, ccode, 3);
    GetAndPrintConst("epshqr2", "epshqr2", 1, 1.0E-20, ccode, 3);
    GetAndPrintConst("epsbisecr", "epsbisecr", 1, 1.0e-10, ccode, 3);
    GetAndPrintConst("epsmach", "epsmach", 1, 1e-15, ccode, 3);
    GetAndPrintConst("epsdigit", "epsdigit", 1, 1e-14, ccode, 3);
    GetAndPrintConst("prnscr", "print_on_screen", 0, 1.0, ccode, 3);
    GetAndPrintConst("prnfile", "print_on_file", 0, 1.0, ccode, 3);
    GetAndPrintConst("failprn", "failure_print", 0, 0.0, ccode, 3);
    GetAndPrintConst("paramregu", "para_regU", 0, 1.0, ccode, 3);
    GetAndPrintConst("qccalc", "qc_calc", 0, 0.0, ccode, 3);
    GetAndPrintConst("hfcalc", "hf_calc", 0, 1.0, ccode, 3);
    GetAndPrintConst("dfcalc", "df_calc", 0, 0.0, ccode, 3);
    GetAndPrintConst("hasscalc", "hass_calc", 0, 1.0, ccode, 3);
    GetAndPrintConst("kuracalc", "kura_calc", 0, 0.0, ccode, 3);
    GetAndPrintConst("vecprn", "eigvec_prn", 0, 0.0, ccode, 3);

    GetStrConst("datafile", name);
    if (GetError()==NotFound)
      fprintf(ccode, "name_datafile= 'kinwrkdat';\n");
    else
      fprintf(ccode, "name_datafile= '%s';\n", name);
    GetStrConst("datafile", name);
    if (GetError()==NotFound)
      fprintf(ccode, "name_textfile= 'kinwrkdat.t';\n");
    else
      fprintf(ccode, "name_textfile= '%s.t';\n", name);

    dyn=0;
    for(i=1; i<=NoOfSpec(); i++) {
      GetSpecNo(i, name, &charge);
      if (IsSpecInConstraint(name, charge)==0) {
	dyn++;
	fprintf(ccode, "xx[%d]=%e;\n", dyn, GetBeginConc(name, charge));

	temp=GetSpecConst(name, charge, "Bfstep");
	BfErrorCode= GetError();
	if (BfErrorCode==NotFound) {
	  fprintf(ccode, "initndir[%d]=-1;\n", dyn);
	  fprintf(ccode, "inithmax[%d]=%e;\n", dyn,1.0E2);
        } else {
	  if (BfErrorCode==NoError) {
	     if (temp>0.0)
	       fprintf(ccode, "initndir[%d]=1;\n", dyn);
	     else
	       fprintf(ccode, "initndir[%d]=-1;\n", dyn);
	     fprintf(ccode, "inithmax[%d]=%e;\n", dyn, fabs(temp));
	  } /* if */
	} /* else */

	temp=GetSpecConst(name, charge, "Bfxmax");
	BfErrorCode= GetError();
	if (BfErrorCode==NotFound) {
	  fprintf(ccode, "initxupp[%d]=%e;\n", dyn, 1.0E3);
        } else {
	  if (BfErrorCode==NoError) {
       	     fprintf(ccode, "initxupp[%d]=%e;\n", dyn, temp);
          } /* if */
        } /* else */

	temp=GetSpecConst(name, charge, "Bfxmin");
	BfErrorCode= GetError();
	if (BfErrorCode==NotFound) {
	  fprintf(ccode, "initxlow[%d]=%e;\n", dyn, 0.0);
        } else {
	  if (BfErrorCode==NoError) {
	     fprintf(ccode, "initxlow[%d]=%e;\n", dyn, temp);
          } /* if */
        } /* else */

	temp=GetSpecConst(name, charge, "Bfpref");
	BfErrorCode= GetError();
	if (BfErrorCode==NotFound) {
	  fprintf(ccode, "initpref[%d]=%e;\n", dyn, 0.1);
        } else {
	  if (BfErrorCode==NoError) {
	     fprintf(ccode, "initpref[%d]=%e;\n", dyn, temp);
          } /* if */
        } /* else */

      } /* if */
    } /* for i */
    for(i=1; i<=NumOfParameter(); i++) {
      GetParamNo(i, name, &charge, &form);
      if (form==1) {
	GetInitParam(name, &temp);
	fprintf(ccode, "xx[n%d]=%e;\n", i, temp);
	temp=GetDirectForParam(name);
	if (temp>0.0)
	  fprintf(ccode, "initndir[n%d]=1;\n", i);
	else
	  fprintf(ccode, "initndir[n%d]=-1;\n", i);
	GetDeltaParam(name, &temp);
	fprintf(ccode, "inithmax[n%d]=%e;\n", i, fabs(temp));
	GetLowHighPrefParam(name, &temp, &temp1, &temp2);
	fprintf(ccode, "initxlow[n%d]=%e;\n", i, temp);
	fprintf(ccode, "initxupp[n%d]=%e;\n", i, temp1);
	if ((temp2>1.0) && (temp2<0.0))
	  fprintf(ccode, "KnCont: Pref. for %s must between 0 and 1.\n", name);
	else
	  fprintf(ccode, "initpref[n%d]=%e;\n", i, temp2);
      } /* if */
      else {
	GetInitConc(name, charge, &temp);
	fprintf(ccode, "xx[n%d]=%e;\n", i, temp);
	temp=GetDirectForConc(name, charge);
	if (temp>0.0)
	  fprintf(ccode, "initndir[n%d]=1;\n", i);
	else
	  fprintf(ccode, "initndir[n%d]=-1;\n", i);
	GetDeltaConc(name, charge, &temp);
	fprintf(ccode, "inithmax[n%d]=%e;\n", i, fabs(temp));
	GetLowHighPrefConc(name, charge, &temp, &temp1, &temp2);
	fprintf(ccode, "initxlow[n%d]=%e;\n", i, temp);
	fprintf(ccode, "initxupp[n%d]=%e;\n", i, temp1);
	if ((temp2>1.0) && (temp2<0.0))
	  fprintf(ccode, "KnCont: Pref. for %s(%d) must between 0 and 1.\n",
		  name, (int)charge);
	else
	  fprintf(ccode, "initpref[n%d]=%e;\n", i, temp2);
      } /* else */
    } /* for i */
    fprintf(ccode, "} /* hopfinit */\n\n");

    /* End of printing hopfinit */

    break;
  } /* switch NumbOfParams */

  /* Printing detnumparam */

  fprintf(ccode, "void detnumparam(void)\n");
  fprintf(ccode, "{\n");
  fprintf(ccode, "  numparam=%d;\n", NumbOfParams);
  fprintf(ccode, "} /* detnumparam */\n\n");
  /* End of printing detnumparam */


  StringFree(name);
  StringFree(rename);
} /* KnCC */
