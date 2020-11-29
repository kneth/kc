/**************************************************************************** 
  Implementation of symbolic mathematic.
  
  (C) Copyright 1992-1996 by
  Kenneth Geisshirt (kneth@fatou.ruc.dk)            Keld Nielsen (kn@kiku.dk)
  Dept. of Life Sciences and Chemistry         Dept. of Theoretical Chemistry
  Roskilde University                                University of Copenhagen
  P.O. Box 5                                             Universitetsparken 5
  4000 Roskilde                                          2100 Copenhagen East
  Denmark                                                             Denmark

  See kc.tex for details.

  Last updated: 17 February 1995
****************************************************************************/

#include "config.h"
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "symbmath.h" 
#include "misc.h"

#ifdef _USE_GARBAGE_COL_
#  include <gc.h>
#else                    /* use ordinary malloc */
#  include <stdlib.h>
#endif

#undef WORKING

/****************************************************************************
  TreeReduce tries to reduce the complexity of the expression. The routine
  is an internal one only, i.e. it is not called by the user explicitely.

  The code is _not_ commented, but at many reductions one can find a short
  explaination. The notation is: Expression 1 -> Expression 2. It should be
  read as: Expression 1 is reduced to expression 2. The symbols f and g 
  represents two generic expressions.
*****************************************************************************/

void TreeReduce(Tree tin) {

  Tree   res, t;
  int    alter;

  alter=0;
  res=TreeCreate();
  t=TreeCreate();
  if (tin!=NULL) {
    switch (tin->tag) {
    case Oper:
      alter=2;
      TreeCpy(t, tin);
      TreeReduce(t->data.oper.left);
      TreeReduce(t->data.oper.right);
      switch (t->data.oper.op) {
      case Add: 
	if ((t->data.oper.left->tag==Konst) && 
	    (t->data.oper.right->tag==Konst)) {
	  alter=1;
	  res->tag=Konst;
	  res->data.value=t->data.oper.left->data.value
	    +t->data.oper.right->data.value;
	} /* if */
	if ((t->data.oper.right->tag==Konst) && (alter==2)) {
	  if (t->data.oper.right->data.value==0.0) {
	    alter=1;
	    TreeCpy(res, t->data.oper.left);
	  } /* if */
	} /* if */
	if ((t->data.oper.left->tag==Konst) && (alter==2)) {
	  if (t->data.oper.left->data.value==0.0) {
	    alter=1;
	    TreeCpy(res, t->data.oper.right);
	  } /* if */
	} /* if */

	/* (f+(-g)) -> (f-g) */
	if ((t->data.oper.right->tag==NegSign) && (alter==2)) {
	  alter=1;
	  res->tag=Oper;
	  res->data.oper.op=Sub;
	  res->data.oper.left=TreeCreate();
	  TreeCpy(res->data.oper.left, t->data.oper.left);
	  res->data.oper.right=TreeCreate();
	  TreeCpy(res->data.oper.right, t->data.oper.right);
	} /* if */
	break;
      case Sub: 
	if ((t->data.oper.left->tag==Konst) && 
	    (t->data.oper.right->tag==Konst)) {
	  alter=1;
	  res->tag=Konst;
	  res->data.value=t->data.oper.left->data.value
	    -t->data.oper.right->data.value;
	} /* if */
	if ((t->data.oper.right->tag==Konst) && (alter==2)) 
	  if (t->data.oper.right->data.value==0.0) {
	    alter=1;
	    TreeCpy(res, t->data.oper.left);
	  } /* if */

	if ((t->data.oper.left->tag==Konst) && (alter==2)) {
	  if (t->data.oper.left->data.value==0.0) {
	    alter=1;
	    res->tag=NegSign;
	    res->data.right=TreeCreate();
	    TreeCpy(res->data.right, t->data.oper.right);
	  }
	} /* if */
	if ((t->data.oper.left->tag==Var) && 
	    (t->data.oper.right->tag==Var)) 
	  if (strcmp(t->data.oper.left->data.name, 
		     t->data.oper.right->data.name)==0) {
	    alter=1;
	    res->tag=Konst;
	    res->data.value=0.0;
	  } /* if */

	/* (f-(-g)) -> (f+g) */
	if ((t->data.oper.right->tag==NegSign) && (alter==2)) {
	  alter=1;
	  res->tag=Oper;
	  res->data.oper.op=Add;
	  res->data.oper.left=TreeCreate();
	  TreeCpy(res->data.oper.left, t->data.oper.left);
	  res->data.oper.right=TreeCreate();
	  TreeCpy(res->data.oper.right, t->data.oper.right);
	} /* if */
	
	/* (-f-g) -> -(f+g) */
	if ((t->data.oper.left->tag==NegSign) && (alter==2)) {
	  alter=1;
	  res->tag=NegSign;
	  res->data.right=TreeCreate();
	  res->data.right->tag=Oper;
	  res->data.right->data.oper.op=Add;
	  res->data.right->data.oper.left=TreeCreate();
	  TreeCpy(res->data.right->data.oper.left, 
		  t->data.oper.left->data.right);
	  res->data.right->data.oper.right=TreeCreate();
	  TreeCpy(res->data.right->data.oper.right, t->data.oper.right);
	}
	break;

      case Mul: 

	/* (Konst1*Konst2) -> Konst3 */
	if ((t->data.oper.left->tag==Konst) && 
	    (t->data.oper.right->tag==Konst)) {
	  alter=1;
	  res->tag=Konst;
	  res->data.value=t->data.oper.left->data.value
	    *t->data.oper.right->data.value;
	} /* if */

	/* (f*0) -> 0 and (f*1) -> (f) */
	if ((t->data.oper.right->tag==Konst) && (alter==2)) {
	  if (t->data.oper.right->data.value==0.0) {
	    alter=1;
	    res->tag=Konst;
	    res->data.value=0.0;
	  }; /* if */
	  if (t->data.oper.right->data.value==1.0) {
	    alter=1;
	    TreeCpy(res, t->data.oper.left);
	  } /* if */
	} /* if */

	/* (0*f) -> 0 and (1*f) -> (f) */
	if ((t->data.oper.left->tag==Konst) && (alter==2)) {
	  if (t->data.oper.left->data.value==0.0) {
	    alter=1;
	    res->tag=Konst;
	    res->data.value=0.0;
	  } /* if */
	  if (t->data.oper.left->data.value==1.0) {
	    alter=1;
	    TreeCpy(res, t->data.oper.right);
	  } /* if */
	} /* if */

#ifdef WORKING
	/* (-1*f) -> (-f) */
	if ((t->data.oper.left->tag==Konst) && (alter==2)) 
	  if (t->data.oper.left->data.value==-1.0) {
	    alter=1;
	    res->tag=NegSign;
	    res->data.right=TreeCreate();
	    TreeCpy(res->data.right, t->data.oper.right);
	  } /* if */

	/* (f*(-1)) -> (-f) */
	if ((t->data.oper.right->tag==Konst) && (alter==2)) 
	  if (t->data.oper.right->data.value==-1.0) {
	    alter=1;
	    res->tag=NegSign;
	    res->data.right=TreeCreate();
	    TreeCpy(res->data.right, t->data.oper.left);
	  } /* if */
#endif

	break;
      case Div: 
	if ((t->data.oper.left->tag==Konst) && 
	    (t->data.oper.right->tag==Konst)) {
	  /* division by zero is not detected */
	  alter=1;
	  res->tag=Konst;
	  res->data.value=t->data.oper.left->data.value
	    /t->data.oper.right->data.value;
	} /* if */
	if ((t->data.oper.right->tag==Konst) && (alter==2)) {
	  if (t->data.oper.right->data.value==1.0) {
	    alter=1;
	    TreeCpy(res, t->data.oper.left);
	  } /* if */
	  if (t->data.oper.right->data.value==0.0) {
	    alter=1;
	    res->tag=Konst;
	    res->data.value=0.0;
	  } /* if */
	} /* if */
	if ((t->data.oper.left->tag==Var) && (t->data.oper.right->tag==Var)) {
	  if (strcmp(t->data.oper.left->data.name, 
		     t->data.oper.right->data.name)==0) {
	    alter=1;
	    res->tag=Konst;
	    res->data.value=1.0;
	  } /* if */
	} /* if */
	if ((t->data.oper.left->tag==Konst) && (alter==2))
	  if (t->data.oper.left->data.value==0.0) {
	    alter=1;
	    res->tag=Konst;
	    res->data.value=0.0;
	  } /* if */
	break;
      case Pow: 
	if ((t->data.oper.left->tag==Konst) && 
	    (t->data.oper.right->tag==Konst)) {
	  alter=1;
	  res->tag=Konst;
	  res->data.value=pow(t->data.oper.left->data.value, 
			      t->data.oper.right->data.value);
	} /* if */
	if ((t->data.oper.right->tag==Konst) && (alter==2)) {
	  if (t->data.oper.right->data.value==1.0) {
	    alter=1;
	    TreeCpy(res, t->data.oper.left);
	  } /* if */
	  if (t->data.oper.right->data.value==0.0) {
	    alter=1;
	    res->tag=Konst;
	    res->data.value=1.0;
	  } /* if */
	} /* if */
	break;
      } /* switch */
      break;
    case NegSign:
      alter=2;
      t=TreeCreate();
      TreeCpy(t, tin);
      TreeReduce(t->data.right);

      /* (-(-f)) -> (f) */
      if (t->data.right->tag==NegSign) {
	alter=1;
	res->tag=t->data.right->data.right->tag;
	res->data.right=TreeCreate();
	TreeCpy(res, t->data.right->data.right);
      } /* if */

      /* -(Constant) -> (-Constant) */
      if (t->data.right->tag==Konst) {
	alter=1;
	res->tag=Konst;
	res->data.value=-t->data.right->data.value;
      } /* if */
      break;
    case Func:
      alter=2;
      t=TreeCreate();
      TreeCpy(t, tin);
      TreeReduce(t->data.func.right);
      break;
    default:
      /* Nothing */
      alter=0;
      break;
    }; /* switch tag */
    switch (alter) {
    case 0: /* nothing altered */
      TreeCpy(res, tin);
      break;
    case 1: /* res is the result */ 
      break;
    case 2: /* t is the result */
      TreeCpy(res, t); 
      break;
    }; /* switch */
/*
    TreeKill(tin);
    tin=TreeCreate();
*/
    TreeCpy(tin, res);
    TreeKill(res);
    TreeKill(t);
  } /* if */
  else {
    fprintf(stderr, "ERROR in TreeReduce: NULL pointer.\nAborting kc - sorry.\n");
    exit(-1);
  }
} /* TreeReduce */


/****************************************************************************
  The following two procdures are used by TreePrint.
*****************************************************************************/

void SplitFP(double x, double *mantissa, int *exponent) {

  int i, sign;
  char str[STRING_LENGTH];

  *mantissa=0.0;
  sprintf(str, "%e", x);
  if (str[0]=='-')
    sign=-1;
  else
    sign=1;
  i=1;
  while (str[i]!='e') { 
    *mantissa=*mantissa+(double)str[i]*pow(10.0, -(double)(i-1));
    i++;
  } /* while */
  *mantissa=*mantissa*sign;
  i+=2;
  if (str[i]=='-') 
    sign=-1;
  else
    sign=1;
  *exponent=10*(int)str[i+1];
  *exponent+=(int)str[i+2];
  *exponent*=sign;
} /* SplitFP */


void LineBreak(FILE *out, int mode) {

  if ((len>=60)) {
    fprintf(out, "\n");
    switch (mode) {
      case 1: 
        fprintf(out, "     &");
        break;
      case 2:
        break;
      case 3:
        /* nothing */
        break;
      } /* switch mode */
    len=0;
  } /* if */
} /* LineBreak */


/****************************************************************************
  The following routines are exported.
*************************************************************************** */


int TreeGetError() {

  return tree_error;
} /* TreeGetError */


Tree TreeCreate(void) {

  Tree t;

#ifdef _USE_GARBAGE_COL_
  t=(Tree) GC_malloc(sizeof(TreeNode));
#else
  t=(Tree) malloc(sizeof(TreeNode));
#endif
  if (t==NULL) {
    fprintf(stderr, "ERROR in TreeCreate: NULL pointer recieved.\nAborting kc - sorry.\n");
    exit(-1);
  }
  tree_error=0;
  return t;
} /* TreeCreate */


void TreeAdd(Tree t1, Tree t2) {
  
  Tree    t;

  t=TreeCreate();
  t->tag=Oper;
  t->data.oper.op=Add;
  t->data.oper.left=TreeCreate();
  t->data.oper.right=TreeCreate();
  TreeCpy(t->data.oper.left, t1);
  TreeCpy(t->data.oper.right, t2);
  tree_error=0;
  TreeReduce(t);
  TreeCpy(t1, t);
  TreeKill(t);
} /* TreeAdd */

void TreeSub(Tree t1, Tree t2) {
  
  Tree    t;

  t=TreeCreate();
  t->tag=Oper;
  t->data.oper.op=Sub;
  t->data.oper.left=TreeCreate();
  t->data.oper.right=TreeCreate();
  TreeCpy(t->data.oper.left, t1);
  TreeCpy(t->data.oper.right, t2);
  tree_error=0;
  TreeReduce(t);
  TreeCpy(t1, t);
  TreeKill(t);
} /* TreeSub */ 

void TreeMul(Tree t1, Tree t2) {
  
  Tree    t;

  t=TreeCreate();
  t->tag=Oper;
  t->data.oper.op=Mul;
  t->data.oper.left=TreeCreate();
  t->data.oper.right=TreeCreate();
  TreeCpy(t->data.oper.left, t1);
  TreeCpy(t->data.oper.right, t2);
  tree_error=0;
  TreeReduce(t);
  TreeCpy(t1, t);
  TreeKill(t);
} /* TreeMul */ 

void TreeDiv(Tree t1, Tree t2) {
  
  Tree    t;

  t=TreeCreate();
  t->tag=Oper;
  t->data.oper.op=Div;
  t->data.oper.left=TreeCreate();
  t->data.oper.right=TreeCreate();
  TreeCpy(t->data.oper.left, t1);
  TreeCpy(t->data.oper.right, t2);
  tree_error=0;
  TreeReduce(t); 
  TreeCpy(t1, t);
  TreeKill(t);
} /* TreeDiv */

void TreePow(Tree t, Tree t2) {

  Tree     res;

  res=TreeCreate();
  res->tag=Oper;
  res->data.oper.op=Pow;
  res->data.oper.left=TreeCreate();
  res->data.oper.right=TreeCreate();
  TreeCpy(res->data.oper.left, t);
  TreeCpy(res->data.oper.right, t2);
  TreeReduce(res);
  tree_error=0;
  TreeCpy(t, res);
  TreeKill(res);
} /* TreePow */

void TreeSign(Tree t) {

  Tree res;

  res=TreeCreate();
  res->tag=NegSign;
  res->data.right=TreeCreate();
  TreeCpy(res->data.right, t);
  TreeReduce(res);
/*  TreeKill(t); 
  t=TreeCreate(); */
  TreeCpy(t, res);
  TreeKill(res);
  tree_error=0;
} /* TreeSign */

void TreeAssignConst(Tree t, double val) {

  t->tag=Konst;
  t->data.value=val;
  tree_error=0;
} /* TreeAssignConst */


void TreeAssignVar(Tree t, char *name) {

  t->tag=Var;
  (void) strcpy(t->data.name, name);
  tree_error=0;
} /* TreeAssignVar */


void TreeSubstVar(Tree t, char* name, double val) {

  if (t==NULL) 
    tree_error=1; 
  else { 
    switch (t->tag) {
    case Var:  
      if (strcmp(t->data.name, name)==0) {
	t->tag=Konst;
	t->data.value=val;
	tree_error=0;
      } /* if */
      break;
    case Oper:
      TreeSubstVar(t->data.oper.left, name, val);
      TreeSubstVar(t->data.oper.right, name, val);
      break;
    case NegSign:
      TreeSubstVar(t->data.right, name, val);
      break;
    case Func:
      TreeSubstVar(t->data.func.right, name, val);
      break;
    } /* switch tag */
    tree_error=0;
  } /* else */
  TreeReduce(t);
} /* TreeSubstVar */

void TreeSubstTree(Tree t, char *name, Tree value) {

  if (t==NULL) 
    tree_error=1; 
  else { 
    switch (t->tag) {
    case Var:
      if (strcmp(t->data.name, name)==0)
	TreeCpy(t, value);
      break;
    case Oper:
      TreeSubstTree(t->data.oper.left, name, value);
      TreeSubstTree(t->data.oper.right, name, value);
      break;
    case Func:
      TreeSubstTree(t->data.func.right, name, value);
      break;
    case NegSign:
      TreeSubstTree(t->data.right, name, value);
      break;
    } /* switch */
    tree_error=0;
    TreeReduce(t);
  } /* else */
} /* TreeSubstTree */

void TreeDerive(Tree res, Tree t, char* name) {

  Tree t1, t2, t3;

  tree_error=0;
  if (t!=NULL) {
    switch (t->tag) {
    case Konst: 
      res->tag=Konst;
      res->data.value=0.0; 
      break;
    case Var: 
      res->tag=Konst;
      if (strcmp(t->data.name, name)==0)
	res->data.value=1.0;
      else
	res->data.value=0.0;
      break;
    case Oper:
      switch (t->data.oper.op) {
      case Add:
	res->tag=Oper;
        res->data.oper.op=Add;
        res->data.oper.left=TreeCreate();
        res->data.oper.right=TreeCreate();
        TreeDerive(res->data.oper.left, t->data.oper.left, name);
        TreeDerive(res->data.oper.right, t->data.oper.right, name);
        break;
      case Sub: 
	res->tag=Oper;           
	res->data.oper.op=Sub;
	res->data.oper.left=TreeCreate();
	res->data.oper.right=TreeCreate();
	TreeDerive(res->data.oper.left, t->data.oper.left, name);
	TreeDerive(res->data.oper.right, t->data.oper.right, name);
	break;
      case Mul: 
	res->data.oper.left=TreeCreate();
	res->data.oper.right=TreeCreate();
	res->tag=Oper;
	res->data.oper.op=Add;
	res->data.oper.left->tag=Oper;
	res->data.oper.left->data.oper.op=Mul;
	res->data.oper.right->tag=Oper;
	res->data.oper.right->data.oper.op=Mul;
	res->data.oper.right->data.oper.left=TreeCreate();
	res->data.oper.right->data.oper.right=TreeCreate();
	res->data.oper.left->data.oper.left=TreeCreate();
	res->data.oper.left->data.oper.right=TreeCreate();
	TreeDerive(res->data.oper.left->data.oper.left, 
		   t->data.oper.left, name);
	TreeCpy(res->data.oper.left->data.oper.right, t->data.oper.right);
	TreeCpy(res->data.oper.right->data.oper.left, t->data.oper.left);
	TreeDerive(res->data.oper.right->data.oper.right, 
		   t->data.oper.right, name);
	break;
      case Div: 
	res->tag=Oper;
	res->data.oper.op=Div;
	res->data.oper.left=TreeCreate();
	res->data.oper.right=TreeCreate();
	res->data.oper.right->tag=Oper;
	res->data.oper.right->data.oper.op=Pow;
	res->data.oper.right->data.oper.left=TreeCreate();
	res->data.oper.right->data.oper.right=TreeCreate();
	res->data.oper.right->data.oper.right->tag=Konst;
	res->data.oper.right->data.oper.right->data.value=2.0;
	TreeCpy(res->data.oper.right->data.oper.left, t->data.oper.right); 
	res->data.oper.left->tag=Oper;
	res->data.oper.left->data.oper.op=Sub;
	res->data.oper.left->data.oper.left=TreeCreate();
	res->data.oper.left->data.oper.right=TreeCreate();
	res->data.oper.left->data.oper.left->tag=Oper;
	res->data.oper.left->data.oper.left->data.oper.op=Mul;
	res->data.oper.left->data.oper.left->data.oper.left=TreeCreate();
	res->data.oper.left->data.oper.left->data.oper.right=TreeCreate();
	TreeCpy(res->data.oper.left->data.oper.left->data.oper.left, 
		t->data.oper.right);
	TreeDerive(res->data.oper.left->data.oper.left->data.oper.right, 
		   t->data.oper.left, name);
	res->data.oper.left->data.oper.right->tag=Oper;
	res->data.oper.left->data.oper.right->data.oper.op=Mul;
	res->data.oper.left->data.oper.right->data.oper.left=TreeCreate();
	res->data.oper.left->data.oper.right->data.oper.right=TreeCreate();
	TreeCpy(res->data.oper.left->data.oper.right->data.oper.left, 
		t->data.oper.left);
	TreeDerive(res->data.oper.left->data.oper.right->data.oper.right, 
		   t->data.oper.right, name);
	break;
      case Pow:             /* (f^g)' = g*f^(g-1)*f' + f^g*log(f)*g' */
	t1=TreeCreate();
	t2=TreeCreate();
	TreeAssignConst(t1, 1.0);
	TreeCpy(t2, t->data.oper.right);
	TreeSub(t2, t1);
	TreeKill(t1);
	t1=TreeCreate();
	TreeCpy(t1, t->data.oper.left);
	TreePow(t1, t2);
	TreeKill(t2);
	t2=TreeCreate();
	TreeDerive(t2, t->data.oper.left, name);
	TreeMul(t1, t2);
	TreeMul(t1, t->data.oper.right);       /* t1 is now first term */
	TreeKill(t2);
	t2=TreeCreate();
	
	TreeCpy(t2, t->data.oper.left);
	TreePow(t2, t->data.oper.right);
	t3=TreeCreate();
	TreeCpy(t3, t->data.oper.left);
	TreeApplyFunc(&t3, Ln);
	TreeMul(t2, t3);
	TreeKill(t3);
	t3=TreeCreate();
	TreeDerive(t3, t->data.oper.right, name);
	TreeMul(t2, t3);                     /* t2 is now second term */
	TreeKill(t3);       
	TreeCpy(res, t1);
	TreeKill(t1);
	TreeAdd(res, t2);                    /* res is the final result */
	TreeKill(t2);
	break;
      } /* switch */
      break;
    case Func:
      res->tag=Oper;
      res->data.oper.op=Mul;
      res->data.oper.left=TreeCreate();
      res->data.oper.right=TreeCreate();
      TreeDerive(res->data.oper.left, t->data.func.right, name);
      switch (t->data.func.fun) {
      case Exp:
	res->data.oper.right->tag=Func;
	res->data.oper.right->data.func.right=TreeCreate();
	TreeCpy(res->data.oper.right->data.func.right, t->data.func.right);
	res->data.oper.right->data.func.fun=Exp;
	break;
      case Cos:
	res->data.oper.right->tag=NegSign;
	res->data.oper.right->data.right=TreeCreate();
	res->data.oper.right->data.right->tag=Func;
	res->data.oper.right->data.right->data.func.right=TreeCreate();
	res->data.oper.right->data.right->data.func.fun=Sin;
	TreeCpy(res->data.oper.right->data.right->data.func.right, 
		t->data.func.right);
	break;
      case Sin:
	res->data.oper.right->tag=Func;
	res->data.oper.right->data.func.fun=Cos;
	res->data.oper.right->data.func.right=TreeCreate();
	TreeCpy(res->data.oper.right->data.func.right, t->data.func.right); 
	break;
      case Cosh:
	res->data.oper.right->tag=Func;
	res->data.oper.right->data.func.fun=Sinh;
	res->data.oper.right->data.func.right=TreeCreate();
	TreeCpy(res->data.oper.right->data.func.right, t->data.func.right);
	break;
      case Sinh:
        res->data.oper.right->tag=Func;
        res->data.oper.right->data.func.fun=Cosh;
        res->data.oper.right->data.func.right=TreeCreate();
        TreeCpy(res->data.oper.right->data.func.right, t->data.func.right);
        break;
      case Atan:
	res->data.oper.op=Div;
	TreeDerive(res->data.oper.left, t->data.oper.right, name);
	res->data.oper.right->tag=Oper;
	res->data.oper.right->data.oper.op=Add;
	res->data.oper.right->data.oper.left=TreeCreate();
	res->data.oper.right->data.oper.right=TreeCreate();
	res->data.oper.right->data.oper.left->tag=Konst;
	res->data.oper.right->data.oper.left->data.value=1.0;
	res->data.oper.right->data.oper.right->tag=Oper;
	res->data.oper.right->data.oper.right->data.oper.op=Pow;
	res->data.oper.right->data.oper.right->data.oper.left=TreeCreate();
	res->data.oper.right->data.oper.right->data.oper.right=TreeCreate();
	TreeCpy(res->data.oper.right->data.oper.right->data.oper.left, 
		t->data.func.right);
	res->data.oper.right->data.oper.right->data.oper.right->tag=Konst;
	res->data.oper.right->data.oper.right->data.oper.right->data.value
	  =2.0;
	break;
      case Asinh:
	res->data.oper.op=Div;
	TreeDerive(res->data.oper.left, t->data.oper.right, name);
	res->data.oper.right->tag=Oper;
	res->data.oper.right->data.oper.op=Pow;
	res->data.oper.right->data.oper.left=TreeCreate();
	res->data.oper.right->data.oper.right=TreeCreate();
	res->data.oper.right->data.oper.left->tag=Oper;
	res->data.oper.right->data.oper.left->data.oper.op=Pow;
	res->data.oper.right->data.oper.left->data.oper.left=TreeCreate();
	res->data.oper.right->data.oper.left->data.oper.right=TreeCreate();
	res->data.oper.right->data.oper.left->data.oper.left->tag=Oper;
	res->data.oper.right->data.oper.left->data.oper.left->data.oper.op
	  =Add;
	res->data.oper.right->data.oper.left->data.oper.left->data.oper.left
	  =TreeCreate();
	res->data.oper.right->data.oper.left->data.oper.left->data.oper.right
	  =TreeCreate();
	res->data.oper.right->data.oper.left->data.oper.left->data.oper.left->tag=Oper;
	res->data.oper.right->data.oper.left->data.oper.left->data.oper.left->data.oper.op=Pow;
	res->data.oper.right->data.oper.left->data.oper.left->data.oper.left->data.oper.left=TreeCreate();
	res->data.oper.right->data.oper.left->data.oper.left->data.oper.left->data.oper.right=TreeCreate();
	TreeCpy(res->data.oper.right->data.oper.left->data.oper.left->data.oper.left->data.oper.left, t->data.func.right);
	res->data.oper.right->data.oper.left->data.oper.left->data.oper.left->data.oper.right->tag=Konst;
	res->data.oper.right->data.oper.left->data.oper.left->data.oper.left->data.oper.right->data.value=2.0;
	res->data.oper.right->data.oper.left->data.oper.left->data.oper.right->tag=Konst;
	res->data.oper.right->data.oper.left->data.oper.left->data.oper.right->data.value=1.0;
	res->data.oper.right->data.oper.left->data.oper.right->tag=Konst;
	res->data.oper.right->data.oper.left->data.oper.right->data.value=0.5;
	break;
      case Atanh:
        res->data.oper.op=Div;
        TreeDerive(res->data.oper.left, t->data.oper.right, name);
        res->data.oper.right->tag=Oper;
        res->data.oper.right->data.oper.op=Sub;
        res->data.oper.right->data.oper.left=TreeCreate();
        res->data.oper.right->data.oper.right=TreeCreate();
        res->data.oper.right->data.oper.left->tag=Konst;
        res->data.oper.right->data.oper.left->data.value=1.0;
        res->data.oper.right->data.oper.right->tag=Oper;
        res->data.oper.right->data.oper.right->data.oper.op=Pow;
        res->data.oper.right->data.oper.right->data.oper.left=TreeCreate();
        res->data.oper.right->data.oper.right->data.oper.right=TreeCreate();
        TreeCpy(res->data.oper.right->data.oper.right->data.oper.left, 
		t->data.func.right);
        res->data.oper.right->data.oper.right->data.oper.right->tag=Konst;
        res->data.oper.right->data.oper.right->data.oper.right->data.value=2.0;
        break;
      case Tan:
	res->data.oper.op=Div;
        TreeDerive(res->data.oper.left, t->data.oper.right, name);
	res->data.oper.right->tag=Oper;
	res->data.oper.right->data.oper.op=Pow;
	res->data.oper.right->data.oper.right=TreeCreate();
	res->data.oper.right->data.oper.left=TreeCreate();
	res->data.oper.right->data.oper.left->tag=Func;
	res->data.oper.right->data.func.fun=Cos;
	res->data.oper.right->data.func.right=TreeCreate();
	TreeCpy(res->data.oper.right->data.func.right, t->data.func.right);
	res->data.oper.right->data.oper.right->tag=Konst;
	res->data.oper.right->data.oper.right->data.value=2.0;
	break;
      case Tanh:
        res->data.oper.op=Div;
        TreeDerive(res->data.oper.left, t->data.oper.right, name);
        res->data.oper.right->tag=Oper;
        res->data.oper.right->data.oper.op=Pow;
        res->data.oper.right->data.oper.right=TreeCreate();
        res->data.oper.right->data.oper.left=TreeCreate();
        res->data.oper.right->data.oper.left->tag=Func;
        res->data.oper.right->data.func.fun=Cosh;
        res->data.oper.right->data.func.right=TreeCreate();
        TreeCpy(res->data.oper.right->data.func.right, t->data.func.right);
        res->data.oper.right->data.oper.right->tag=Konst;
        res->data.oper.right->data.oper.right->data.value=2.0;
	break;
      case Asin:
	res->data.oper.op=Div;
	res->data.oper.left=TreeCreate();
	res->data.oper.right=TreeCreate();
	TreeDerive(res->data.oper.left, t->data.func.right, name);
	res->data.oper.right->tag=Oper;
	res->data.oper.right->data.oper.op=Pow;
	res->data.oper.right->data.oper.left=TreeCreate();
	res->data.oper.right->data.oper.right=TreeCreate();
	res->data.oper.right->data.oper.left->tag=Oper;
	res->data.oper.right->data.oper.left->data.oper.op=Sub;
	res->data.oper.right->data.oper.left->data.oper.left=TreeCreate();
	res->data.oper.right->data.oper.left->data.oper.right=TreeCreate();
	res->data.oper.right->data.oper.left->data.oper.left->tag=Konst;
	res->data.oper.right->data.oper.left->data.oper.left->data.value=1.0;
	res->data.oper.right->data.oper.left->data.oper.right->tag=Oper;
	res->data.oper.right->data.oper.left->data.oper.right->data.oper.op=Pow;
	res->data.oper.right->data.oper.left->data.oper.right->data.oper.left=TreeCreate();
	res->data.oper.right->data.oper.left->data.oper.right->data.oper.right=TreeCreate();
	TreeCpy(res->data.oper.right->data.oper.left->data.oper.right->data.oper.left, t->data.func.right);
	res->data.oper.right->data.oper.left->data.oper.right->data.oper.right->tag=Konst;
	res->data.oper.right->data.oper.left->data.oper.right->data.oper.right->data.value=2.0;
	res->data.oper.right->data.oper.right->tag=Konst;
	res->data.oper.right->data.oper.right->data.value=0.5;
	break;
      case Acos:
	res->data.oper.op=Div;
	res->data.oper.left=TreeCreate();
	res->data.oper.right=TreeCreate();
	res->data.oper.left->tag=NegSign;
	TreeDerive(res->data.oper.left->data.right, t->data.func.right, name);
	res->data.oper.right->tag=Oper;
	res->data.oper.right->data.oper.op=Pow;
	res->data.oper.right->data.oper.left=TreeCreate();
	res->data.oper.right->data.oper.right=TreeCreate();
	res->data.oper.right->data.oper.left->tag=Oper;
	res->data.oper.right->data.oper.left->data.oper.op=Sub;
	res->data.oper.right->data.oper.left->data.oper.left=TreeCreate();
	res->data.oper.right->data.oper.left->data.oper.right=TreeCreate();
	res->data.oper.right->data.oper.left->data.oper.left->tag=Konst;
	res->data.oper.right->data.oper.left->data.oper.left->data.value=1.0;
	res->data.oper.right->data.oper.left->data.oper.right->tag=Oper;
	res->data.oper.right->data.oper.left->data.oper.right->data.oper.op=Pow;
	res->data.oper.right->data.oper.left->data.oper.right->data.oper.left=TreeCreate();
	res->data.oper.right->data.oper.left->data.oper.right->data.oper.right=TreeCreate();
	TreeCpy(res->data.oper.right->data.oper.left->data.oper.right->data.oper.left, t->data.func.right);
	res->data.oper.right->data.oper.left->data.oper.right->data.oper.right->tag=Konst;
	res->data.oper.right->data.oper.left->data.oper.right->data.oper.right->data.value=2.0;
	res->data.oper.right->data.oper.right->tag=Konst;
	res->data.oper.right->data.oper.right->data.value=0.5;
	break;
      case Acosh:
	fprintf(stderr, "TreeDerive: dAcosh/dx not implemented\n");
	tree_error=4;
	break;
      case Log:
	res->data.oper.right->tag=Oper;
	res->data.oper.right->data.oper.op=Div;
	res->data.oper.right->data.oper.left=TreeCreate();
	res->data.oper.right->data.oper.right=TreeCreate();
	res->data.oper.right->data.oper.left->tag=Konst;
        res->data.oper.right->data.oper.left->data.value=0.4342944819;
        TreeCpy(res->data.oper.right->data.oper.right, t->data.func.right);
	break;
      }; /* switch */
      break; /* Functions */
    case NegSign:
      res->tag=NegSign;
      res->data.right=TreeCreate();
      TreeDerive(res->data.right, t->data.right, name);
      break;
    }; /* switch */
    TreeReduce(res); 
  }; /* if */
} /* TreeDerive */

double TreeEval(Tree t) {
 
  double c1, c2, res;

  if (t==NULL) 
    return 0.0;
  else {
    switch (t->tag) {
    case Var:
      tree_error=NoEval;
      res=0.0;
      break;
    case Konst:
      res=t->data.value;
      tree_error=TreeNoError;
      break;
    case Oper: 
      c1=TreeEval(t->data.oper.left);
      if (tree_error==TreeNoError) {
	c2=TreeEval(t->data.oper.right);
	if (tree_error==TreeNoError) {
	  switch (t->data.oper.op) {
	  case Add: 
	    res=c1+c2; 
	    break;
	  case Sub: 
	    res=c1-c2; 
	    break;
	  case Mul: 
	    res=c1*c2; 
	    break;
	  case Div: 
	    res=c1/c2; 
	    break;
	  case Pow: 
	    res=pow(c1, c2); 
	    break;
	  }; /* switch */
	  tree_error=TreeNoError;
	} /* if */
      } /* if */
      break;
    case NegSign:
      res=-TreeEval(t->data.right);
      break;
    case Func:
      res=TreeEval(t->data.func.right);
      if (tree_error==TreeNoError) {
	switch (t->data.func.fun) {
	case Exp:
	  res=exp(res);
	  break;
	case Ln:
	  res=log(res);
	  break;
	case Log:
	  res=log10(res);
	  break;
	case Sin:
	  res=sin(res);
	  break;
	case Cos:
	  res=cos(res);
	  break;
	case Tan:
	  res=tan(res);
	  break;
	case Sinh:
	  res=sinh(res);
	  break;
	case Cosh:
	  res=cosh(res);
	  break;
	case Tanh:
	  res=tanh(res);
	  break;
	case Asin:
	  res=asin(res);
	  break;
	case Acos:
	  res=acos(res);
	  break;
	case Atan:
	  res=atan(res);
	  break;
	case Asinh:
	  /* res=asinh(res); */
	  break;
	case Acosh:
	  /* res=acosh(res); */
	  break;
	case Atanh:
	  /* res=atanh(res); */
	  break;
	} /* switch */
      } /* if */
      break;
    } /* switch */
  } /* else */
 return res;
} /* TreeEval */

/* The following procedure is 'internal' */

void TreePrintXtra(Tree t, int mode, FILE *output) {

/* Mode:  1 - (HP) Fortran 77 */
/*	  2 - (HP) Pascal     */
/*        3 - ANSI C          */

  int i, expo;
  double mant;

  switch (mode) {
  case 1: 
    switch (t->tag) {
    case Oper:  
      fprintf(output, "("); 
      len++;
      LineBreak(output, mode);
      TreePrintXtra(t->data.oper.left, mode, output);
      switch (t->data.oper.op) {
      case Add: 
	fprintf(output, "+"); 
	len++; 
	LineBreak(output, mode); 
	break;
      case Sub: 
	fprintf(output, "-"); 
	len++; 
	LineBreak(output, mode); 
	break;
      case Div: 
	fprintf(output, "/"); 
	len++; 
	LineBreak(output, mode); 
	break;
      case Mul: 
	fprintf(output, "*"); 
	len++; 
	LineBreak(output, mode); 
	break;
      case Pow: 
	fprintf(output, "**"); 
	len=len+2; 
	LineBreak(output, mode); 
	break;
      }; /* switch */
      TreePrintXtra(t->data.oper.right, mode, output);
      fprintf(output, ")"); 
      len++;
      LineBreak(output, mode);
      break;
    case Konst:
      fprintf(output, "(%e)", t->data.value);
      len+=15;
      /*
      if (fmod(t->data.value, 1.0)==0.0) {
	fprintf(output, "(%d)", (int)t->data.value);
	len+=4;
      } else {
	fprintf(output, "(%e)", t->data.value);
	len+=15;
      }*/ /* if */
      LineBreak(output, mode);
      break;
    case Var:   
      len=len+strlen(t->data.name); 
      fprintf(output, "%s", t->data.name); 
      LineBreak(output, mode);
      break; 
    case Func:
      switch (t->data.func.fun) {
      case Exp:
	fprintf(output, "EXP(");
	len+=4;
	break;
      case Log:
	fprintf(output, "LOG10(");
	len+=6;
	break;
      case Ln:
	fprintf(output, "LOG(");
	len+=4;
	break;
      case Sin:
	fprintf(output, "SIN(");
	len+=4;
	break;
      case Cos:
	fprintf(output, "COS(");
	len+=4;
	break;
      case Asinh:
	fprintf(output, "ASINH(");
	len+=6;
	break;
      case Acosh:
	fprintf(output, "ACOSH(");
	len+=6;
	break;
      case Atanh:
	fprintf(output, "ATANH(");
	len+=6;
	break;
      case Tan:
	fprintf(output, "TAN(");
	len+=4;
	break;
      } /* switch */
      LineBreak(output, mode);
      TreePrintXtra(t->data.func.right, mode, output);
      fprintf(output, ")");
      len++;
      LineBreak(output, mode);
      break; /* Func */
    case NegSign:
      fprintf(output, "(-(");
      len+=2;
      LineBreak(output, mode);
      TreePrintXtra(t->data.right, mode, output);
      fprintf(output, "))");
      len++;
      LineBreak(output, mode);
      break; /* NegSign */
    }; /* switch */
    break;
  case 2: 
    switch (t->tag) {
    case Oper:  
      fprintf(output, "("); 
      len++;
      LineBreak(output, mode);
      switch (t->data.oper.op) {
      case Add: 
	TreePrintXtra(t->data.oper.left, mode, output);
	fprintf(output, "+");
	len++;
	LineBreak(output, mode);
	TreePrintXtra(t->data.oper.right, mode, output);
	break;
      case Sub: 
	TreePrintXtra(t->data.oper.left, mode, output);
	fprintf(output, "-"); 
	len++;
	LineBreak(output, mode);
	TreePrintXtra(t->data.oper.right, mode, output);
	break;
      case Div: 
	TreePrintXtra(t->data.oper.left, mode, output);
	fprintf(output, "/"); 
	len++;
	LineBreak(output, mode);
	TreePrintXtra(t->data.oper.right, mode, output);
	break;
      case Mul: 
	TreePrintXtra(t->data.oper.left, mode, output);
	fprintf(output, "*"); 
	len++;
	LineBreak(output, mode);
	TreePrintXtra(t->data.oper.right, mode, output);
	break;
      case Pow: 
	if (t->data.oper.right->tag==Konst) {
	  if (ceil(t->data.oper.right->data.value)==floor(t->data.oper.right->data.value)) {
	    TreePrintXtra(t->data.oper.left, mode, output);
	    for(i=1; i<((int)(t->data.oper.right->data.value)); i++) {
	      fprintf(output, "*");
	      TreePrintXtra(t->data.oper.left, mode, output);
	      len++;
	      LineBreak(output, mode);
	    }
	  } else {
	    fprintf(output, "EXP("); /* Pascal has no power operator!! */ 
	    len+=4;
	    LineBreak(output, mode);
	    TreePrintXtra(t->data.oper.right, mode, output);
	    fprintf(output, "*LOG(");
	    len+=5;
	    LineBreak(output, mode);
	    TreePrintXtra(t->data.oper.left, mode, output);
	    fprintf(output, "))");
	    len+=2;
	    LineBreak(output, mode);
	  } /* if ... else ... */
	} else {
	  fprintf(output, "EXP("); /* Pascal has no power operator!! */ 
	  len+=4;
	  LineBreak(output, mode);
	  TreePrintXtra(t->data.oper.right, mode, output);
	  fprintf(output, "*LOG(");
	  len+=5;
	  LineBreak(output, mode);
	  TreePrintXtra(t->data.oper.left, mode, output);
	  fprintf(output, "))");
	  len+=2;
	  LineBreak(output, mode);
	}
	break;
      } /* switch */
      fprintf(output, ")"); 
      len++;
      LineBreak(output, mode);
      break;
    case Konst:
      if (fmod(t->data.value, 1.0)==0.0) {
	fprintf(output, "(%d)", (int)t->data.value);
	len+=5;
      } else {
	fprintf(output, "(%e)", t->data.value);
	len+=15;
      } /* if */
      LineBreak(output, mode);
      break;
    case Var:    
      fprintf(output, "%s", t->data.name); 
      len+=strlen(t->data.name);
      LineBreak(output, mode);
      break;
    case Func:
      switch (t->data.func.fun) {
      case Exp:
        fprintf(output, "exp(");
        len+=4;
        break;
      case Log:
        fprintf(output, "log10(");
        len+=6;
        break;
      case Ln:
        fprintf(output, "log(");
        len+=4;
        break;
      case Sin:
        fprintf(output, "sin(");
        len+=4;
        break;
      case Cos:
        fprintf(output, "cos(");
        len+=4;
        break;
      case Tan:
        fprintf(output, "tan(");
        len+=4;
        break;
      case Sinh:
        fprintf(output, "sinh(");
        len+=5;
        break;
      case Cosh:
        fprintf(output, "cosh(");
        len+=5;
        break;
      case Tanh:
        fprintf(output, "tanh(");
        len+=5;
        break;
      case Asin:
        fprintf(output, "asin(");
        len+=5;
        break;
      case Acos:
        fprintf(output, "acos(");
        len+=5;
        break;
      case Atan:
        fprintf(output, "atan(");
        len+=5;
        break;
      case Asinh:
        fprintf(output, "asinh(");
        len+=6;
        break;
      case Acosh:
        fprintf(output, "acosh(");
        len+=6;
        break;
      case Atanh:
        fprintf(output, "atanh(");
        len+=6;
        break;
      } /* switch */
      LineBreak(output, mode);
      TreePrintXtra(t->data.func.right, mode, output);
      fprintf(output, ")");
      len++;
      LineBreak(output, mode);
      break; /* Func */
    case NegSign:
      fprintf(output, "(-(");
      len+=2;
      LineBreak(output, mode);
      TreePrintXtra(t->data.right, mode, output);
      fprintf(output, "))");
      len++;
      LineBreak(output, mode);
      break; /* NegSign */
    }; /* switch */
    break; 
  case 3: /* The C mode */
    switch (t->tag) {
    case Var:   
      fprintf(output, "%s", t->data.name);
      len+=strlen(t->data.name);
      LineBreak(output, mode);
      break;
    case Konst: 
      fprintf(output, "%e", t->data.value);
      len+=15;
      LineBreak(output, mode);
      break;
    case Oper:  
      fprintf(output, "(");
      len++;
      LineBreak(output, mode);
      switch (t->data.oper.op) {
      case Div: 
	TreePrintXtra(t->data.oper.left, mode, output);
	fprintf(output, "/");
	len++;
	LineBreak(output, mode);
	TreePrintXtra(t->data.oper.right, mode, output);
	break;
      case Mul: 
	TreePrintXtra(t->data.oper.left, mode, output);
	fprintf(output, "*");
	len++;
	LineBreak(output, mode);
	TreePrintXtra(t->data.oper.right, mode, output);
	break;
      case Add: 
	TreePrintXtra(t->data.oper.left, mode, output);
	fprintf(output, "+");
	len++;
	LineBreak(output, mode);
	TreePrintXtra(t->data.oper.right, mode, output);
	break;
      case Sub:
	TreePrintXtra(t->data.oper.left, mode, output);
	fprintf(output, "-");
	len++;
	LineBreak(output, mode);
	TreePrintXtra(t->data.oper.right, mode, output);
	break;
      case Pow:
	if (t->data.oper.right->tag==Konst) {
	  if (ceil(t->data.oper.right->data.value)==floor(t->data.oper.right->data.value)) {
	    TreePrintXtra(t->data.oper.left, mode, output);
	    for(i=1; i<((int)(ceil(t->data.oper.right->data.value))); i++) {
	      fprintf(output, "*");
	      len++;
	      LineBreak(output, mode);
	      TreePrintXtra(t->data.oper.left, mode, output);
	    } /* for i */
	  } else {
	    fprintf(output, "pow(");
	    len+=4;
	    LineBreak(output, mode);
	    TreePrintXtra(t->data.oper.left, mode, output);
	    fprintf(output, ",");
	    len++;
	    LineBreak(output, mode);
	    TreePrintXtra(t->data.oper.right, mode, output);
	    fprintf(output, ")");
	    len++;
	    LineBreak(output, mode);
	  } /* if ... else ... */
	} else {
	  fprintf(output, "pow(");
	  len+=4;
	  LineBreak(output, mode);
	  TreePrintXtra(t->data.oper.left, mode, output);
	  fprintf(output, ",");
	  len++;
	  LineBreak(output, mode);
	  TreePrintXtra(t->data.oper.right, mode, output);
	  fprintf(output, ")");
	  len++;
	  LineBreak(output, mode);
	} /* if ... else ... */
	break;
      } /* switch oper */
      fprintf(output, ")");
      len++;
      LineBreak(output, mode);
      break; /* Oper */
    case Func:
      switch (t->data.func.fun) {
      case Exp:
	fprintf(output, "exp(");
	len+=4;
	break;
      case Log:
	fprintf(output, "log10(");
	len+=6;
	break;
      case Ln:
	fprintf(output, "log(");
	len+=4;
	break;
      case Sin:
	fprintf(output, "sin(");
	len+=4;
	break;
      case Cos:
	fprintf(output, "cos(");
	len+=4;
	break;
      case Tan:
	fprintf(output, "tan(");
	len+=4;
	break;
      case Sinh:
	fprintf(output, "sinh(");
	len+=5;
	break;
      case Cosh:
	fprintf(output, "cosh(");
	len+=5;
	break;
      case Tanh:
	fprintf(output, "tanh(");
	len+=5;
	break;
      case Asin:
	fprintf(output, "asin(");
	len+=5;
	break;
      case Acos:
	fprintf(output, "acos(");
	len+=5;
	break;
      case Atan:
	fprintf(output, "atan(");
	len+=5;
	break;
      case Asinh:
	fprintf(output, "asinh(");
	len+=6;
	break;
      case Acosh:
	fprintf(output, "acosh(");
	len+=6;
	break;
      case Atanh:
	fprintf(output, "atanh(");
	len+=6;
	break;
      } /* switch */
      LineBreak(output, mode);
      TreePrintXtra(t->data.func.right, mode, output);
      fprintf(output, ")");
      len++;
      LineBreak(output, mode);
      break; /* Func */
    case NegSign:
      fprintf(output, "(-(");
      len+=2;
      LineBreak(output, mode);
      TreePrintXtra(t->data.right, mode, output);
      fprintf(output, "))");
      len++;
      LineBreak(output, mode);
      break; /* NegSign */
    }; /* switch tag */
    break; /* mode 3 */   
  }; /* switch */
} /* TreePrintXtra */


void TreePrint(Tree t, int mode, FILE *output) {

  len=10;
  TreePrintXtra(t, mode, output);
} /* TreePrint */

void TreeCpy(Tree res, Tree t) {  
  
  if (t==NULL) {
    res=NULL;
    tree_error=NoTree;
  } /* if */
  else {
    switch (t->tag) {
    case Oper:  
      res->tag=Oper;
      res->data.oper.left=TreeCreate();
      res->data.oper.right=TreeCreate();
      res->data.oper.op=t->data.oper.op;
      TreeCpy(res->data.oper.left, t->data.oper.left);
      TreeCpy(res->data.oper.right, t->data.oper.right);
      break;
    case Konst: 
      res->tag=Konst;
      res->data.value=t->data.value;
      break;
    case Var:   
      res->tag=Var;
      (void) strcpy(res->data.name, t->data.name);
      break;
    case Func:
      res->tag=Func;
      res->data.func.fun=t->data.func.fun;
      res->data.func.right=TreeCreate();
      TreeCpy(res->data.func.right, t->data.func.right);
      break;
    case NegSign:
      res->tag=NegSign;
      res->data.right=TreeCreate();
      TreeCpy(res->data.right, t->data.right);
      break;
     }; /* switch */
  }; /* else */
} /* TreeCpy */

void TreeKill(Tree t) {

#ifndef _USE_GARBAGE_COL_
  if (t!=NULL) { 
    switch (t->tag) {
    case Oper: 
      TreeKill(t->data.oper.left);
      TreeKill(t->data.oper.right);
      break; 
    case Var: 
      break;
    case Func:
      TreeKill(t->data.func.right);
      break;
    case NegSign:
      TreeKill(t->data.right);
      break;
    } 
    free((MALLOCTYPE *) t);
  }  
#endif
} /* TreeKill */

void TreeApplyFunc(Tree *t, Function fun) {

  Tree temp;

  temp=TreeCreate();
  temp->tag=Func;
  temp->data.func.fun=fun;
  temp->data.func.right=TreeCreate();
  TreeCpy(temp->data.func.right, *t);
  TreeReduce(temp);
  TreeKill(*t);
  *t=TreeCreate(); 
  TreeCpy(*t, temp);
  TreeKill(temp);
  tree_error=0;
} /* TreeApplyFunc */  
