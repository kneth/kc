/************************************************************************** 
  CharMode - a code generator for kc and the CharMode
  package. 

  CopyWrong 1993-1994 by
  Kenneth Gesshirt (kneth@osc.kiku.dk)
  Department of Theoretical Chemistry
  H.C. Orsted Institute
  Universitetsparken 5
  2100 Copenhagen
  Denmark

  See kc.tex for details.
  Last updated: 12 May 1995  by KN
***************************************************************************/

#include <stdio.h>
#include <malloc.h>
#include <time.h>
#include <string.h>
#include "config.h"
#include "symbmath.h"
#include "tableman.h"
#include "codegen.h"
#include "misc.h"

void CharMode(FILE* ccode, FILE *hcode) {

  double  charge, temp, coeff;
  char    *name, *rename, *strtmp;
  time_t  timer;
  Tree    v_temp, tmp, temp_tree, tree_temp;
  int     i, j, NumbOfDynVars, react_no, finished, constraint, dyn, dyn2;

  name=StringAlloc();
  rename=StringAlloc();

  NumbOfDynVars= NoOfSpec()+NumOfDynVar()-NumOfConstraint();

  InitCodeGenVar(NumbOfDynVars, NumOfConstraint(), NoOfReact());
  GenerateRateExpr();
  GenerateJacobi();
  GenerateDiffReac();

  timer=time(&timer);
  fprintf(hcode, "/*********************************************\n");
  fprintf(hcode, " WARNING: This file was generated by kc v%s\n", VERSION);
  fprintf(hcode, " CopyWrong by Kenneth Geisshirt\n");
  fprintf(hcode, " %s", ctime(&timer));
  fprintf(hcode, "**********************************************/\n");
  fprintf(hcode, "#ifndef _MODEL_HEADER_\n#define _MODEL_HEADER_\n");
  fprintf(hcode, "#define novar_ %d\n", NumbOfDynVars);
  fprintf(hcode, "#define noreac_ %d\n", NoOfReact());
  fprintf(hcode, "double x_[novar_];\n");
  fprintf(hcode, "int do_print[novar_];\n");
  fprintf(hcode, "char species[25][novar_];\n");
  fprintf(hcode, "double rfw_[noreac_];\n");
  fprintf(hcode, "double rrv_[noreac_];\n");
  fprintf(hcode, "double v_[novar_];\n");
  fprintf(hcode, "double rfwds_[noreac_][novar_];\n");
  fprintf(hcode, "double rrvds_[noreac_][novar_];\n");
  fprintf(hcode, "double jacobi_matx[novar_][novar_];\n");
  fprintf(hcode, "extern void ReacRate(double *);\n");
  fprintf(hcode, "extern void ReacFlow(double *);\n");
  fprintf(hcode, "extern void Jacobi(double *);\n");
  fprintf(hcode, "extern void ReacFlow_ds(double *);\n");
  fprintf(hcode, "extern void InitValues(void);\n");
  fprintf(hcode, "#endif\n");



  fprintf(ccode, "/*********************************************\n");
  fprintf(ccode, " WARNING: This file was generated by kc v%s\n", VERSION);
  fprintf(ccode, " CopyWrong by Kenneth Geisshirt\n");
  fprintf(ccode, " %s", ctime(&timer));
  fprintf(ccode, "**********************************************/\n");
  fprintf(ccode, "#include \"model.h\"\n");
  fprintf(ccode, "#include <math.h>\n\n");

  /* Printing ReacRate */

  fprintf(ccode, "void ReacRate(double* S_) {\n");
    for(i=1; i<=NoOfSpec(); i++) {
      GetSpecNo(i, name, &charge);
      if (IsSpecInConstraint(name, charge)==0) {
	RenameSpec(rename, name, charge);
	fprintf(ccode, "double %s;\n", rename);
      } /* if */
    } /* for i */
    for(i=1; i<=NumOfDynVar(); i++) {
      GetDynVarNo(i, name);
      fprintf(ccode, "double %s;\n", name);
    } /* for i */
    dyn=1;
    for(i=1; i<=NoOfSpec(); i++) {
      GetSpecNo(i, name, &charge);
      if (IsSpecInConstraint(name, charge)==0) {
	RenameSpec(rename, name, charge);
	fprintf(ccode, "%s=S_[%d];\n", rename, dyn-1);
	dyn++;
      } /* if */
    } /* for i*/
    for(i=1; i<=NumOfDynVar(); i++) {
      GetDynVarNo(i, name);
      fprintf(ccode, "%s=S_[%d];\n", name, i+NoOfSpec()-NumOfConstraint()-1);
    } /* for i */
  dyn=0;
  for(i=1; i<=(NoOfSpec()+NumOfDynVar()-NumOfConstraint()); i++) {
    fprintf(ccode, "v_[%d]=(", i-1);
    TreePrint(v[i-1], 3, ccode);
    fprintf(ccode, ");\n");
  } /* for i */
  fprintf(ccode, "} /* ReacRate */\n\n");

  /* End of printing ReacRate */

  /* Printing Jacobi matrix routine */

  fprintf(ccode, "void Jacobi(double *S_){\n");
    for(i=1; i<=NoOfSpec(); i++) {
      GetSpecNo(i, name, &charge);
      if (IsSpecInConstraint(name, charge)==0) {
	RenameSpec(rename, name, charge);
	fprintf(ccode, "double %s;\n", rename);
      } /* if */
    } /* for i */
    for(i=1; i<=NumOfDynVar(); i++) {
      GetDynVarNo(i, name);
      fprintf(ccode, "double %s;\n", name);
    } /* for i */
    dyn=1;
    for(i=1; i<=NoOfSpec(); i++) {
      GetSpecNo(i, name, &charge);
      if (IsSpecInConstraint(name, charge)==0) {
	RenameSpec(rename, name, charge);
	fprintf(ccode, "%s=S_[%d];\n", rename, dyn-1);
	dyn++;
      } /* if */
    } /* for i*/
    for(i=1; i<=NumOfDynVar(); i++) {
      GetDynVarNo(i, name);
      fprintf(ccode, "%s=S_[%d];\n", name, i+NoOfSpec()-NumOfConstraint()-1);
    } /* for i */

  for(i=0; i<NumbOfDynVars; i++)
    for(j=0; j<NumbOfDynVars; j++) {
      temp=TreeEval(jacobi[i][j]);
      if (TreeGetError()==NoEval) {
	fprintf(ccode, "jacobi_matx[%d][%d]=", i, j);
	TreePrint(jacobi[i][j], 3, ccode);
	fprintf(ccode, ";\n");
      } /* if */
    } /* for j */
  fprintf(ccode, "} /* Jacobi */\n\n");

  /* End of printing jacobimatrix */

  /* Printing ReacFlow */

  fprintf(ccode, "void ReacFlow(double* S_) {\n");
    for(i=1; i<=NoOfSpec(); i++) {
      GetSpecNo(i, name, &charge);
      if (IsSpecInConstraint(name, charge)==0) {
	RenameSpec(rename, name, charge);
	fprintf(ccode, "double %s;\n", rename);
      } /* if */
    } /* for i */
    for(i=1; i<=NumOfDynVar(); i++) {
      GetDynVarNo(i, name);
      fprintf(ccode, "double %s;\n", name);
    } /* for i */
    dyn=1;
    for(i=1; i<=NoOfSpec(); i++) {
      GetSpecNo(i, name, &charge);
      if (IsSpecInConstraint(name, charge)==0) {
	RenameSpec(rename, name, charge);
	fprintf(ccode, "%s=S_[%d];\n", rename, dyn-1);
	dyn++;
      } /* if */
    } /* for i*/
    for(i=1; i<=NumOfDynVar(); i++) {
      GetDynVarNo(i, name);
      fprintf(ccode, "%s=S_[%d];\n", name, i+NoOfSpec()-NumOfConstraint()-1);
    } /* for i */
  for(i=1; i<=NoOfReact(); i++) {
    fprintf(ccode, "rfw_[%d]=(", i-1);
    TreePrint(rfw[i-1], 3, ccode);
    fprintf(ccode, ");\n");
    fprintf(ccode, "rrv_[%d]=(", i-1);
    TreePrint(rrv[i-1], 3, ccode);
    fprintf(ccode, ");\n");
  } /* for i */
  fprintf(ccode, "} /* ReacFlow */\n\n");

  /* End of printing ReacRate */

  /* Printing ReacFlowdS */

  fprintf(ccode, "void ReacFlowdS(double* S_) {\n");
    for(i=1; i<=NoOfSpec(); i++) {
      GetSpecNo(i, name, &charge);
      if (IsSpecInConstraint(name, charge)==0) {
	RenameSpec(rename, name, charge);
	fprintf(ccode, "double %s;\n", rename);
      } /* if */
    } /* for i */
    for(i=1; i<=NumOfDynVar(); i++) {
      GetDynVarNo(i, name);
      fprintf(ccode, "double %s;\n", name);
    } /* for i */
    dyn=1;
    for(i=1; i<=NoOfSpec(); i++) {
      GetSpecNo(i, name, &charge);
      if (IsSpecInConstraint(name, charge)==0) {
	RenameSpec(rename, name, charge);
	fprintf(ccode, "%s=S_[%d];\n", rename, dyn-1);
	dyn++;
      } /* if */
    } /* for i*/
    for(i=1; i<=NumOfDynVar(); i++) {
      GetDynVarNo(i, name);
      fprintf(ccode, "%s=S_[%d];\n", name, i+NoOfSpec()-NumOfConstraint()-1);
    } /* for i */
  for(i=1; i<=NoOfReact(); i++) {
    for(j=1; j<=NumbOfDynVars; j++) {
       fprintf(ccode, "rfwds_[%d][%d]=(", i-1,j-1);
       TreePrint(rfwds[i-1][j-1], 3, ccode);
       fprintf(ccode, ");\n");
       fprintf(ccode, "rrvds_[%d][%d]=(", i-1,j-1);
       TreePrint(rrvds[i-1][j-1], 3, ccode);
       fprintf(ccode, ");\n");
    }; /* for j */
  } /* for i */
  fprintf(ccode, "} /* ReacFlowdS */\n\n");

  /* End of printing ReacFlowdS */

  /* print init routine */

  fprintf(ccode, "void InitValues(void) {\n");
  dyn=0;
  for(i=0; i<NoOfSpec(); i++) {
    GetSpecNo(i+1, name, &charge);
    if (IsSpecInConstraint(name, charge)==0) {
      temp=GetBeginConc(name, charge);
      if (GetError()==NoError) 
	  fprintf(ccode, "x_[%d]=%e;\n", dyn, temp);
	else
	  fprintf(ccode, "x_[%d]=0.0;\n", dyn);
	dyn++;
      } /* if */
    } /* for i */
    for(i=1; i<=NumOfDynVar(); i++) {
      GetDynVarNo(i, name);
      temp=GetInitValue(name);
      fprintf(ccode, "x_[%d]=%e;\n", i+NoOfSpec()-NumOfConstraint()-1, temp);
    } /* for i */
  GetAndPrintConst("epsr", "epsr_", 1, 1e-5, ccode, 3);
  GetAndPrintConst("epsa", "epsa_", 1, 1e-10, ccode, 3);
  GetAndPrintConst("ref", "ref_", 0, 1, ccode, 3);

    dyn=0;   
    for(i=1; i<=NoOfSpec(); i++) {
      GetSpecNo(i, name, &charge);
      if (IsSpecInConstraint(name, charge)==0) {
	RenameSpec(rename, name, charge);
	fprintf(ccode, "(void) strcpy(species[%d], \"%s\");\n", dyn, rename);
	dyn++;
      } /* if */
    } /* for i */
    for(i=1; i<=NumOfDynVar(); i++) {
      GetDynVarNo(i, name);
      fprintf(ccode, "(void) strcpy(species[%d], \"%s\");\n", 
	      i+NoOfSpec()-NumOfConstraint()-IsNonAutoSystem()-1, name);
    } /* for i */
    if (IsNonAutoSystem()==1)
      fprintf(ccode, "x_[equa-1]=stime_;\n");
    
  GetStrConst("datafile", name);
  if (GetError()==NotFound) 
    fprintf(ccode, "(void) strcpy(datafilename_,\"kinwrkda\");\n");
  else
    fprintf(ccode, "(void) strcpy(datafilename_,\"%s\");\n", name);

  /* constant elements of the jacobi */
    for(i=0; i<NoOfSpec()+NumOfDynVar()-NumOfConstraint(); i++)
      for(j=0; j<NoOfSpec()+NumOfDynVar()-NumOfConstraint(); j++) {
	temp=TreeEval(jacobi[i][j]);
	if (TreeGetError()==TreeNoError) 
	  fprintf(ccode, "jacobi_matx[%d][%d]=%e;\n", i, j, temp);
      } /* for j */  
  fprintf(ccode, "}\n");

  /* End of printing InitValues */

  StringFree(name);
  StringFree(rename);
} /* CharMode */

