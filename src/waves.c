/************************************************************************* 
  Waves - a code generator for kc and KGadi/per1d. 
  This is a full operating generator.   

  CopyWrong 1993-1996 by 
  Kenneth Geisshirt (kneth@fatou.ruc.dk)     Keld Nielsen (kn@kin.kiku.dk)
  Dept. of Life Sciences and Chemistry      Dept. of Theoretical Chemistry
  Roskilde University                             University of Copenhagen
  P.O. Box 260                                        Universitetsparken 5
  4000 Roskilde                                            2100 Copenhagen
  Denmark                                                          Denmark

  Last updated: 17 April 1996 by KG
**************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "config.h"
#include "symbmath.h"
#include "tableman.h"
#include "codegen.h"
#include "misc.h"

void AddFullstops(char* str) {

/* Add fullstops (.) so length(str)=16 */

  int i;
  for(i=strlen(str); i<16; i++) 
    str=strcat(str, ".");
} /* AddFullstops */

void Waves(FILE* code_h, FILE* code_c, FILE* code_ini) {

  double  charge, temp, coeff;
  char    name[STRING_LENGTH], rename[STRING_LENGTH];
  time_t  timer;
  Tree    v_temp, tmp, temp_tree;
  int     i, j, react_no, finished, constraint, dyn, dyn2;

  timer=time(&timer);

  fprintf(code_h, "/* %s*/\n\n", ctime(&timer));
  fprintf(code_c, "/* %s*/\n\n", ctime(&timer));
  fprintf(code_c, "#include <math.h>\n");
  fprintf(code_h, "#define n_grids %d\n", (int)GetConstant("ngrid"));
  if (GetError()==NotFound)
    fprintf(stderr, "ERROR: The constant ngrid is not defined!\n");
  fprintf(code_h, "#define m_grids %d\n", (int)GetConstant("mgrid"));
  if (GetError()==NotFound)
    fprintf(stderr, "ERROR: The constant mgrid is not defined!\n");
  i=NoOfSpec()+NumOfDynVar()-NumOfConstraint(); /* abuse of i */
  fprintf(code_h, "#define equa %d\n\n", i);
  fprintf(code_h, "extern void init_diff_const(void);\n");
  fprintf(code_h, "double D[equa];\n");
  fprintf(code_c, "#include \"model.h\"\n");
  fprintf(code_c, "void init_diff_const(void) {\n");
  dyn=1;
  for(i=0; i<NoOfSpec(); i++) {
    GetSpecNo(i+1, name, &charge);
    if (IsSpecInConstraint(name, charge)==0) {
      temp=GetSpecConst(name, charge, "D");
      if (GetError()==NotFound) 
	fprintf(stderr, "ERROR: Diffusion coeff. for %s(%e) not found\n", 
		name, charge);
      else
	fprintf(code_c, "  D[%d] = %e;\n", dyn-1, temp);
      dyn++;
    } /* if */
  } /* for i */
  fprintf(code_c, "} /* init_diff_const */\n\n");
  fprintf(code_h, "extern double reac(double*, int);\n");
  InitCodeGenVar(NoOfSpec()+NumOfDynVar()-NumOfConstraint(),
		 NumOfConstraint(),NoOfReact());
  GenerateRateExpr();
  GenerateJacobi();
  fprintf(code_c, "double reac(double* S, int l) {\n");
  for(i=1; i<=NoOfSpec(); i++) {
    GetSpecNo(i, name, &charge);
    RenameSpec(rename, name, charge);
    fprintf(code_c, " double %s;\n", rename);
  }; /* for i */
  for(i=1; i<=NumOfDynVar(); i++) {
    GetDynVarNo(i, name);
    fprintf(code_c, " double %s;\n", name);
  }; /* for i */
  dyn=1;
  for(i=1; i<=NoOfSpec(); i++) {
    GetSpecNo(i, name, &charge);
    RenameSpec(rename, name, charge);
    constraint=IsSpecInConstraint(name, charge);
    if (constraint==0) {
      fprintf(code_c, "  %s=S[%d];\n", rename, dyn-1);
      dyn++;
    }; /* if */
  }; /* for i*/
  for(i=1; i<=NumOfDynVar(); i++) {
    GetDynVarNo(i, name);
    fprintf(code_c, "  %s=S[%d];\n", name, i+NoOfSpec()-NumOfConstraint());
  }; /* for i */
  for(i=1; i<=NoOfSpec(); i++) {
    GetSpecNo(i, name, &charge);
    RenameSpec(rename, name, charge);
    constraint=IsSpecInConstraint(name, charge);
    if (constraint>0) {  
      fprintf(code_c, "  %s=", rename);
      TreePrint(con[constraint-1], 3, code_c);
      fprintf(code_c, ";\n");
    }; /* if */
  }; /* for i */
  dyn=0;
  fprintf(code_c, "  switch (l) {\n");
  for(i=1; i<=(NoOfSpec()+NumOfDynVar()-NumOfConstraint()); i++) {
    fprintf(code_c, "  case %d: return (", i-1);
    TreePrint(v[i-1], 3, code_c);
    fprintf(code_c, ");\nbreak;\n");
  }; /* for i */
  fprintf(code_c, "  } /* switch */\n");
  fprintf(code_c, "} /* reac */\n");
  
  /* function eval */
  fprintf(code_h, "extern void eval(double *, double *);\n");
  fprintf(code_c, "void eval(double S[equa], double F[equa]) {\n");
  
  for(i=1; i<=NoOfSpec(); i++) {
    GetSpecNo(i, name, &charge);
    RenameSpec(rename, name, charge);
    fprintf(code_c, " double %s;\n", rename);
  }; /* for i */
  for(i=1; i<=NumOfDynVar(); i++) {
    GetDynVarNo(i, name);
    fprintf(code_c, " double %s;\n", name);
  }; /* for i */
  
  dyn=1;
  for(i=1; i<=NoOfSpec(); i++) {
    GetSpecNo(i, name, &charge);
    RenameSpec(rename, name, charge);
    constraint=IsSpecInConstraint(name, charge);
    if (constraint==0) {
      fprintf(code_c, "  if (S[%d]<0.0) S[%d]=0.0;\n", dyn-1, dyn-1);
      fprintf(code_c, "  %s=S[%d];\n", rename, dyn-1);
      dyn++;
    }; /* if */
  }; /* for i*/
  for(i=1; i<=NumOfDynVar(); i++) {
    GetDynVarNo(i, name);
    fprintf(code_c, "  %s=S[%d];\n", name, i+NoOfSpec()-1);
  }; /* for i */

  for(i=1; i<=NoOfSpec(); i++) {
    GetSpecNo(i, name, &charge);
    RenameSpec(rename, name, charge);
    constraint=IsSpecInConstraint(name, charge);
    if (constraint>0) {
      fprintf(code_c, "  %s=", rename);
      TreePrint(con[constraint-1], 3, code_c);
      fprintf(code_c, ";\n");
    }; /* if */
  }; /* for i */

  dyn=0;
  for(i=1; i<=(NoOfSpec()+NumOfDynVar()-NumOfConstraint()); i++) {
    fprintf(code_c, "  F[%d] = ", i-1);
    TreePrint(v[i-1], 3, code_c);
    fprintf(code_c, ";\n");
  } /* for i */
  fprintf(code_c, "} /* eval */\n");

  /* Printing Jacoby matrix */
  fprintf(code_h, "void calc_jac(double S[equa], double tk, double R[equa], double Jac[equa+1][equa+1]);\n");
  fprintf(code_c, "void calc_jac(double S[equa], double tk, double R[equa], double Jac[equa+1][equa+1]) \n{\n");
  for(i=1; i<=NoOfSpec(); i++) {
    GetSpecNo(i, name, &charge);
    RenameSpec(rename, name, charge);
    fprintf(code_c, " double %s;\n", rename);
  }; /* for i */
  for(i=1; i<=NumOfDynVar(); i++) {
    GetDynVarNo(i, name);
    fprintf(code_c, " double %s;\n", name);
  }; /* for i */
  dyn=1;
  for(i=1; i<=NoOfSpec(); i++) {
    GetSpecNo(i, name, &charge);
    RenameSpec(rename, name, charge);
    constraint=IsSpecInConstraint(name, charge);
    if (constraint==0) {
      fprintf(code_c, "  if (S[%d]<0.0) S[%d]=0.0;\n", dyn-1, dyn-1);
      fprintf(code_c, "  %s=S[%d];\n", rename, dyn-1);
      dyn++;
    }; /* if */
  }; /* for i*/
  for(i=1; i<=NumOfDynVar(); i++) {
    GetDynVarNo(i, name);
    fprintf(code_c, "  %s=S[%d];\n", name, i+NoOfSpec()-1);
  }; /* for i */
  for(i=1; i<=NoOfSpec(); i++) {
    GetSpecNo(i, name, &charge);
    RenameSpec(rename, name, charge);
    constraint=IsSpecInConstraint(name, charge);
    if (constraint>0) {
      fprintf(code_c, "  %s=", rename);
      TreePrint(con[constraint-1], 3, code_c);
      fprintf(code_c, ";\n");
    }; /* if */
  }; /* for i */
  for(i=0; i<NoOfSpec()+NumOfDynVar()-NumOfConstraint(); i++)
    for(j=0; j<NoOfSpec()+NumOfDynVar()-NumOfConstraint(); j++) {
      fprintf(code_c, "Jac[%d][%d]=tk*(", i+1, j+1);
      TreePrint(jacobi[i][j], 3, code_c);
      fprintf(code_c, ")");
      if (i==j)
	fprintf(code_c, "-(1+4*R[%d])", i);
      fprintf(code_c, ";\n");
    } /* for j */
  fprintf(code_c, "}\n");

  /* Generate reac_F routine */
  fprintf(code_h, "extern reac_F(double S[equa], double K[equa], double R[equa], double tk, int l);\n");
  fprintf(code_c, "reac_F(double S[equa], double K[equa], double R[equa], double tk, int l) {\n");
  for(i=1; i<=NoOfSpec(); i++) {
    GetSpecNo(i, name, &charge);
    RenameSpec(rename, name, charge);
    fprintf(code_c, " double %s;\n", rename);
  }; /* for i */
  for(i=1; i<=NumOfDynVar(); i++) {
    GetDynVarNo(i, name);
    fprintf(code_c, " double %s;\n", name);
  }; /* for i */
  dyn=1;
  for(i=1; i<=NoOfSpec(); i++) {
    GetSpecNo(i, name, &charge);
    RenameSpec(rename, name, charge);
    constraint=IsSpecInConstraint(name, charge);
    if (constraint==0) {
      fprintf(code_c, "  %s=S[%d];\n", rename, dyn-1);
      dyn++;
    }; /* if */
  }; /* for i*/
  for(i=1; i<=NumOfDynVar(); i++) {
    GetDynVarNo(i, name);
    fprintf(code_c, "  %s=S[%d];\n", name, i+NoOfSpec()-1);
  }; /* for i */
  for(i=1; i<=NoOfSpec(); i++) {
    GetSpecNo(i, name, &charge);
    RenameSpec(rename, name, charge);
    constraint=IsSpecInConstraint(name, charge);
    if (constraint>0) {
      fprintf(code_c, "  %s=", rename);
      TreePrint(con[constraint-1], 3, code_c);
      fprintf(code_c, ";\n");
    }; /* if */
  }; /* for i */
  fprintf(code_c, "switch (l) {\n");
  for(i=1; i<=NoOfSpec()-NumOfConstraint()+NumOfDynVar(); i++) {
    fprintf(code_c, "case %d:\n", i);
    fprintf(code_c, "return (tk*(");
    TreePrint(v[i-1], 3, code_c);
    fprintf(code_c, ")-S[%d]*(1+4*R[%d])+K[%d]);\nbreak;\n", i-1, i-1, i-1);
  } /* for i */
  fprintf(code_c, "}\n}\n");

  /* Printing Jacobi matrix (the right one) */
  fprintf(code_h, "void calc_jac2(double S[equa], double Jac[equa][equa]);\n");
  fprintf(code_c, "void calc_jac2(double S[equa], double Jac[equa][equa]) \n{\n");
  for(i=1; i<=NoOfSpec(); i++) {
    GetSpecNo(i, name, &charge);
    RenameSpec(rename, name, charge);
    fprintf(code_c, "double %s;\n", rename);
  }; /* for i */
  for(i=1; i<=NumOfDynVar(); i++) {
    GetDynVarNo(i, name);
    fprintf(code_c, "double %s;\n", name);
  }; /* for i */
  dyn=1;
  for(i=1; i<=NoOfSpec(); i++) {
    GetSpecNo(i, name, &charge);
    RenameSpec(rename, name, charge);
    constraint=IsSpecInConstraint(name, charge);
    if (constraint==0) {
      fprintf(code_c, "%s=S[%d];\n", rename, dyn-1);
      dyn++;
    }; /* if */
  }; /* for i*/
  for(i=1; i<=NumOfDynVar(); i++) {
    GetDynVarNo(i, name);
    fprintf(code_c, "%s=S[%d];\n", name, i+NoOfSpec()-1);
  }; /* for i */
  for(i=1; i<=NoOfSpec(); i++) {
    GetSpecNo(i, name, &charge);
    RenameSpec(rename, name, charge);
    RenameSpec(rename, name, charge);
    constraint=IsSpecInConstraint(name, charge);
    if (constraint>0) {
      fprintf(code_c, "  %s=", rename);
      TreePrint(con[constraint-1], 3, code_c);
      fprintf(code_c, ";\n");
    }; /* if */
  }; /* for i */
  for(i=0; i<NoOfSpec()+NumOfDynVar()-NumOfConstraint(); i++)
    for(j=0; j<NoOfSpec()+NumOfDynVar()-NumOfConstraint(); j++) {
      fprintf(code_c, "Jac[%d][%d]=", i, j);
      TreePrint(jacobi[i][j], 3, code_c);
      fprintf(code_c, ";\n");
    } /* for j */
  fprintf(code_c, "}\n");

  /* Determine mode */
  temp=GetConstant("mode");
  if (GetError()==NoError)
    switch ((int)temp) {
    case 1: 
      fprintf(code_h, "#define COMPACT\n");
      break;
    case 2:
      fprintf(code_h, "#define DIST\n");
      break;
    case 3:
      fprintf(code_h, "#define ERR_EST\n");
      break;
    case 4:
      fprintf(code_h, "#define SING\n");
      break;
    } /* switch l */
  
  /* Generate input data file */
  temp=GetConstant("dt");
  if (GetError()==NotFound) temp=2.0;
  fprintf(code_ini, "dt......= %e\n", temp);
  temp=GetConstant("length");
  if (GetError()==NotFound) temp=100.0;
  fprintf(code_ini, "L.......= %e\n", temp);
  GetStrConst("prefix", name);
  if (GetError()==NotFound) strcpy(name, "test");
  fprintf(code_ini, "prefix..= %s\n", name);
  temp=GetConstant("tmax");
  if (GetError()==NotFound) temp=200.0;
  fprintf(code_ini, "max_t...= %e\n", temp);
  temp=GetConstant("print1");
  if (GetError()==NotFound) temp=100.0;
  fprintf(code_ini, "pt1.....= %d\n", (int)temp);
  temp=GetConstant("print2");
  if (GetError()==NotFound) temp=100.0;
  fprintf(code_ini, "pt2.....= %d\n", (int)temp);
  temp=GetConstant("update");
  if (GetError()==NotFound) temp=10.0;
  fprintf(code_ini, "period..= %e\n", temp);
  for(i=0; i<NoOfSpec(); i++) {
    GetSpecNo(i, name, &charge);
    if (IsSpecInConstraint(name, charge)==0) {
      temp=GetBeginConc(name, charge);
      fprintf(code_ini, "init-val= %e\n", temp);
    } /* if */
  } /* for i */
  temp=GetConstant("noise");
  if (GetError()==NotFound) temp=0.001;
  fprintf(code_ini, "noise...= %e\n", temp);
} /* Waves */
