%{ 
/********************************************************* 
  Parser to kc.
  CopyWrong 1992-1995 by 
  Kenneth Geisshirt (kneth@fatou.ruc.dk)
  Department of Life Sciences and Chemistry
  Roskilde University
  P.O. Box 260
  4000 Roskilde
  Denmark


  See kc.tex for details.

  Last updated: 17 February 1994
**********************************************************/

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include "symbmath.h"
#include "tableman.h"
#include "codecall.h"

typedef enum {no, ini, ord} Conc; 

typedef struct compound {
  char   name[STRING_LENGTH];
  double charge;
  Conc   concs;
} compound;

typedef struct Strnum {
  int    flag;
  char   name[STRING_LENGTH];
  double numb;
} Strnum;

double     coeff, charge, value, temp, temp1, temp2;
char       flag, name[STRING_LENGTH], string[STRING_LENGTH];
int        i, j, side, lineno=1;
Tree       tmp;
%}

%union {
  double      dval;
  char        oper;
  char        name[STRING_LENGTH];
  compound    comp;
  char        flag;
  Tree        tree;
  Function    func;
  Strnum      strnum;
}

%token <name> names 
%token <name> leftarr
%token <name> rightarr
%token <dval> numbers
%token param 
%token print 
%token <flag> powop
%token semicolon 
%token quotation
%token <flag> equal 
%token <flag> colon 
%token <flag> R 
%token <flag> E 
%token <flag> powc 
%token <flag> leftpar 
%token <flag> rightpar
%token <flag> oneway 
%token <flag> twoways
%token <flag> plus 
%token <flag> minus 
%token <flag> multi 
%token <flag> pdiv
%token <flag> comma 
%token <flag> leftconc 
%token <flag> rightconc 
%token <flag> time0 
%token <flag> K
%token <flag> radical
%token <flag> V
%token <flag> prime 
%token <flag> fun_exp 
%token <flag> fun_log 
%token <flag> fun_ln 
%token <flag> fun_sin 
%token <flag> fun_cos 
%token <flag> fun_tan
%token <flag> fun_sinh 
%token <flag> fun_cosh 
%token <flag> fun_tanh 
%token <flag> fun_asin 
%token <flag> fun_acos 
%token <flag> fun_atan
%token <flag> fun_asinh 
%token <flag> fun_acosh 
%token <flag> fun_atanh
%type <tree> expr 
%type <dval> charge size coeff
%type <oper> sign 
%type <flag> timeopt kind 
%type <comp> subst conc substs
%type <func> function
%type <strnum> strnumb
%start system
%left plus minus
%left multi pdiv
%left powop
%nonassoc UMINUS
%%
system    : { 
	      SetupTableMan(); 
	    }
	    declars reactions constants
	    {
	      if (verbose==1) {               /* Printing status report */
		printf("Species used (*=dyn. var.):\n");
		for(i=1; i<=NoOfSpec(); i++) {
		  GetSpecNo(i, name, &charge);
		  if (IsSpecInConstraint(name, charge)==0)
		    printf("  *%s(", name);
		  else
		    printf("   %s(", name);
		  if (charge==FLT_MAX)
		    printf(".");
		  else {
		    if (charge>0.0)
		      printf("%d+", (int)charge);
		    else if (charge<0.0)
		      printf("%d-", -(int)charge);
		  } /* else */
		  printf(")\n");
		} /* for i*/
		for(i=1; i<=NumOfDynVar(); i++) {
		  GetDynVarNo(i, name);
		  printf("  *%s\n", name);
		} /* for i */
		if (NumOfParameter()!=0) {
		  printf("Parameters declared:\n");
		  for(i=1; i<=NumOfParameter(); i++) {
		    GetParamNo(i, name, &charge, &side); /* abuse of side */
		    if (charge==0.0)
		      printf("  %s\n", name);
		    else {
		      if (side==1) 
			printf("  %s\n", name);
		      else {
			printf("  %s(", name);
			if (charge==FLT_MAX)
			  printf(".");
			else {
			  if (charge>0.0)
			    printf("%d+", (int)charge);
			  else if (charge<0.0)
			    printf("%d-", -(int)charge);
			} /* else */
			printf(")\n");
		      } /* else */
		    } /* else */
		  } /* for i */
		} /* if */
		printf("Constants declared:\n");
		for(i=1; i<=NumOfConstants(); i++) {
		  GetConstantNo(i, name);
		  printf("  %s = %e\n", name, GetConstant(name));
		} /* for i */
		for(i=1; i<=NumOfStrConst(); i++) {
		  GetStrConstNo(i, name, string);
		  printf("  %s = \"%s\"\n", name, string);
		} /* for i */
	      } /* if */
	      if (IsNonAutoSystem()==1) {
		NewDynVar("time");
		tmp=TreeCreate();
		TreeAssignConst(tmp, 1.0);
		NewExpr("time", tmp);
		TreeKill(tmp);
	      } /* if */
	    };
declars   : declars declar
	    | ;
declar    : constant
	    | parameter
	    | print printlist semicolon;
constant  : names equal strnumb 
	    { if ($3.flag==1)
		NewConstant($1, $3.numb);
	      else {
		NewStrConst($1, $3.name);
	      }
	    };
strnumb   : expr semicolon
	    { temp=TreeEval($1);
	      if (TreeGetError()==NoEval) 
		fprintf(stderr, "Unable to evaluate expression in line %d.\n", lineno);
	      else {
		$$.flag=1;
		$$.numb=temp;
	      }
	      TreeKill($1);
	    }
	    | quotation names quotation semicolon
	    { $$.flag=2;
	      strcpy($$.name, $2);
	    };
printlist : printlist comma prn_entry
	    | prn_entry;
prn_entry : names
	    { NewPrintVar($1); }
	    | conc
	    { NewPrintConc($1.name, $1.charge); }; 
parameter : param pentry semicolon;
pentry    : names equal expr comma expr comma expr comma expr comma expr comma expr
	    { temp=TreeEval($3);
	      if (TreeGetError()==NoEval)
		fprintf(stderr, "WARNING: Could not evaluate expr in line %d\n", lineno);
	      else {
		NewParameter($1, temp);
	      }
	      TreeKill($3); 
	      temp=TreeEval($5);
	      if (TreeGetError()==NoEval)
		fprintf(stderr, "WARNING: Could not evaluate expr in line %d\n", lineno);
	      else
		NewDeltaParam($1, temp);
	      TreeKill($5);
	      temp1=TreeEval($7);
	      if (TreeGetError()==NoEval)
		fprintf(stderr, "WARNING: Could not evaluate expr in line %d\n" , lineno);
	      TreeKill($7);
	      temp2=TreeEval($11);
	      if (TreeGetError()==NoEval)
		fprintf(stderr, "WARNING: Could not evaluate expr in line %d\n" , lineno);
	      TreeKill($11);
	      NewLowHighPrefParam($1, temp, temp1, temp2);
	      temp=TreeEval($13);
	      if (TreeGetError()==NoEval)
		fprintf(stderr, "WARNING: Could not evaluate expr in line %d\n", lineno);
	      TreeKill($13);
	      NewDirectForParam($1, (int)temp);
	      temp=TreeEval($9);
	      if (TreeGetError()==NoEval)
		fprintf(stderr, "WARNING: Could not evaluate expr in line %d\n", lineno);
	      TreeKill($9);
	      NewDeltaParam($1, temp);
	    }
	    | conc equal expr comma expr comma expr comma expr comma expr comma expr
	    { temp=TreeEval($3);
	      if (TreeGetError()==NoEval)
		fprintf(stderr, "WARNING: Could not evaluate expr in line %\n", lineno);
	      else
		NewParamConc($1.name, $1.charge, temp);
	      temp=TreeEval($5);
	      if (TreeGetError()==NoEval)
		fprintf(stderr, "WARNING: Could not evaluate expr in line %\n", lineno);
	      TreeKill($3);
	      TreeKill($5);
	      temp1=TreeEval($7);
	      if (TreeGetError()==NoEval)
		fprintf(stderr, "WARNING: Could not evaluate expr in line %d\n" , lineno);
	      TreeKill($7);
	      temp2=TreeEval($11);
	      if (TreeGetError()==NoEval)
		fprintf(stderr, "WARNING: Could not evaluate expr in line %d\n" , lineno);
	      TreeKill($11);
	      NewLowHighPrefConc($1.name, $1.charge, temp, temp1, temp2);
	      temp=TreeEval($13);
	      if (TreeGetError()==NoEval)
		fprintf(stderr, "WARNING: Could not evaluate expr in line %d\n", lineno);
	      TreeKill($9);
	      NewDirectForConc($1.name, $1.charge, (int)temp);
	      temp=TreeEval($9);
	      if (TreeGetError()==NoEval)
		fprintf(stderr, "WARNING: Could not evaluate expr in line %d\n", lineno);
	      TreeKill($9);
	      NewDeltaConc($1.name, $1.charge, temp);
	    };
reactions : reaction
	    | reactions reaction;
reaction  : numbers 
	    { NewReaction((int)value); 
	      side=1; 
	    } 
	    colon substs kind 
	    { 
	      switch ($5) {
	      case '>': 
		AddReactionKind(GetCurrentReaction(), uni); 
		break;
	      case '<': 
		AddReactionKind(GetCurrentReaction(), bi); 
		break;
	      case '=': 
		AddReactionKind(GetCurrentReaction(), equi); 
		break;
	      }; /* switch */
	      side=-1;
	    }
	    substs semicolon reackonst
	    | names prime equal expr semicolon 
	    { NewExpr($1, $4);
	      NewDynVar($1);
	      TreeKill($4);
	    };
reackonst : powconsts semicolon rateconsts  
	    | rateconsts;
powconsts : powconsts semicolon powconst
	    | powconst;
powconst  : names leftpar subst rightpar equal expr
	    { if (strcmp($1, "c")!=0) 
		fprintf(stderr, "c expected!\n");
	      else {
		temp=TreeEval($6);
		if (TreeGetError()==NoEval) 
		  fprintf(stderr, "Unable to evaluate expression in line %d.\n", lineno);
		else
		  NewPowerConst(GetCurrentReaction(), $3.name, $3.charge, temp, side); 
	      }; /* else */
	      TreeKill($6);
	    };
rateconsts: leftarr equal expr  
	    { if (strcmp($1, "k")==0) {  
		NewRateConst(GetCurrentReaction(), -1, $3); 
	      } /* if */
	      else 
		if (strcmp($1, "v")==0) 
		  NewRateExpr(GetCurrentReaction(), -1, $3);
		else
		  fprintf(stderr, "Syntax error: v or k expected in line %d.\n", lineno);
	      TreeKill($3);               
	    } 
	    semicolon ropt
	    | names equal expr 
	    { if (strcmp($1, "K")!=0) 
		fprintf(stderr, "Syntax error: K expected in line %d.\n", lineno);
	      else {
		NewRateConst(GetCurrentReaction(), 0, $3);
	      }; /* else */
	      TreeKill($3); 
	    };
ropt      : { if (GetReactKind(GetCurrentReaction())==bi) (void) fprintf(stderr, "The reaction in line %d is two-ways.\n", lineno);
	    }
	    | rightarr equal expr semicolon
	      { if (GetReactKind(GetCurrentReaction())!=bi)   
		  (void) fprintf(stderr, "The reaction in line %d is a >one-way< reaction or >equilibrium.<\n", lineno);
		else {
		  if (strcmp($1, "k")==0)  
		    NewRateConst(GetCurrentReaction(), 1, $3);  
		  else
		    if (strcmp($1, "v")==0) 
		      NewRateExpr(GetCurrentReaction(), 1, $3);
		    else
		      fprintf(stderr, "Syntax error: v or k expected in line %d.\n", lineno);
	      }; /* else */
	      TreeKill($3);
	    };
kind      : oneway { $$='>'; } 
	    | twoways { $$='<'; }
	    | equal { $$='='; };
substs    : coeff subst    
	    { SpecieInReaction(GetCurrentReaction(), $2.name, $2.charge);
	      NewCoeff(GetCurrentReaction(), $2.name, $2.charge, $1, side);
	      NewSpecie($2.name, $2.charge);
	    }
	    | substs plus coeff subst 
	    { SpecieInReaction(GetCurrentReaction(), $4.name, $4.charge);
	      NewCoeff(GetCurrentReaction(), $4.name, $4.charge, $3, side);
	      NewSpecie($4.name, $4.charge);
	    };
subst     : names charge 
	    { 
	      (void) strcpy($$.name, $1);
	      $$.charge=$2;
	      $$.concs=no;
	    };
charge    : leftpar size rightpar
	    { $$=$2; }
	    | { 
		$$=0.0;
	      }; 
size      : radical 
	    { charge=FLT_MAX;
	      $$=FLT_MAX;
	    }
	    | sign numbers 
	    { if ($1=='+') {
		$$=$2;
		charge=$2;
	      } /* if */
	      else {
		charge=-$2;
		$$=-$2;
	      }; /* else */
	    }
	    | sign
	    { if ($1=='+')
		$$=1.0;
	      else
		$$=-1.0;
	    };
sign      : plus { $$ = '+'; }
	    | minus { $$ = '-'; };
coeff     : numbers
	    { coeff=value;
	      $$=$1;  
	    }
	    | 
	    { coeff=1.0;
	      $$=1.0;
	    };
constants : const  
	    | constants const;
const     : conc timeopt equal expr semicolon
	    { temp=TreeEval($4); 
	      switch ($2) {
	      case 1: 
		if (TreeGetError()==NoEval) 
		  fprintf(stderr, "Unable to evaluate expression.\n");
		else  
		  NewBeginConc($1.name, $1.charge, temp); 
		break;
	      case 0: 
		NewConstraint($1.name, $1.charge, $4);
		break;
	      };
	      TreeKill($4);
	    }
	    | names leftpar subst rightpar equal expr semicolon
	    { temp=TreeEval($6);
	      NewSpecConst($3.name, $3.charge, $1, temp);
	      if (GetError()==NotFound) 
		if ($3.charge==0.0) 
		  NewDynVarConst($3.name, $1, temp);
	      TreeKill($6);
	    };
	    | names timeopt equal expr semicolon
	    { temp=TreeEval($4);
	      if (TreeGetError()==NoEval)
		yyerror("Unable to evaluate expression");
	      else {
		if ($2==1) {
		  NewInitValue($1, temp);
		  TreeKill($4);
		} else
		  yyerror("(0) expected");
	      }
	    }
expr      : minus expr %prec UMINUS
	    { $$=TreeCreate();
	      TreeCpy($$, $2);
	      TreeSign($$);
	      TreeKill($2);
	    }
	    | expr plus expr
	    { $$=TreeCreate(); 
	      TreeAdd($1, $3);
	      TreeCpy($$, $1);
	      TreeKill($3);
	      TreeKill($1);
	    }
	    | expr minus expr
	    { $$=TreeCreate(); 
	      TreeSub($1, $3);
	      TreeCpy($$, $1);
	      TreeKill($3);
	      TreeKill($1);
	    }
	    | expr multi expr
	    {
	      $$=TreeCreate();
	      TreeMul($1, $3);
	      TreeCpy($$, $1);
	      TreeKill($1);
	      TreeKill($3);
	    }
	    | expr pdiv expr
	    {
	      $$=TreeCreate();
	      TreeDiv($1, $3);
	      TreeCpy($$, $1);
	      TreeKill($1);
	      TreeKill($3);
	    }
	    | expr powop expr
	    { $$=TreeCreate();
	      TreeCpy($$, $1);
	      TreePow($$, $3);
	      TreeKill($1);
	      TreeKill($3);
	    }
	    | leftpar expr rightpar 
	    { $$=TreeCreate();
	      TreeCpy($$, $2); 
	      TreeKill($2);
	    }
	    | names 
	    {
	     $$=TreeCreate();
	     temp=GetConstant($1); 
	     if (GetError()==NotFound) {
	       TreeAssignVar($$, $1); 
	       if (strcmp($1, "time")==0)
		 NonAutoSystem();
	     } else
	       TreeAssignConst($$, temp);
	    } 
	    | numbers 
	    { $$=TreeCreate();
	      TreeAssignConst($$, $1);
	    }
	    | conc timeopt
	    {  
	      $$=TreeCreate(); 
	      if ($2==1) {
		temp=GetBeginConc($1.name, $1.charge);
		if (GetError()==NoError)
		  TreeAssignConst($$, temp);
		else
		  fprintf(stderr, "[%s(%e)] not found.\n", $1.name, $1.charge);
		flag='1';
	      } /* if */
	      else { 
		flag='0';
		RenameSpec(name, $1.name, $1.charge);
		TreeAssignVar($$, name);
	      }; /* else */
	    }
	    | function leftpar expr rightpar
	    { $$=TreeCreate();
	      TreeCpy($$, $3);
	      TreeApplyFunc(&$$, $1);
	      TreeKill($3);
	    };
function  : fun_exp     { $$=Exp; }
	    | fun_ln    { $$=Ln; }
	    | fun_log   { $$=Log; }
	    | fun_sin   { $$=Sin; }
	    | fun_cos   { $$=Cos; }
	    | fun_tan   { $$=Tan; }
	    | fun_sinh  { $$=Sinh; }
	    | fun_cosh  { $$=Cosh; }
	    | fun_tanh  { $$=Tanh; }
	    | fun_asin  { $$=Asin; }
	    | fun_acos  { $$=Acos; }
	    | fun_atan  { $$=Atan; }
	    | fun_asinh { $$=Asinh; }
	    | fun_acosh { $$=Acosh; }
	    | fun_atanh { $$=Atanh; }; 
conc      : leftconc subst rightconc
	    { 
	      (void) strcpy($$.name, $2.name);
	      $$.charge=$2.charge;
	      $$.concs=ord;
	    };
timeopt   : { $$=0; } 
	    | time0
	    { $$=1; };
%%
#include "lex.c"
