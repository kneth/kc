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
****************************************************************************/

#include "config.h"
#include <limits.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "symbmath.h" 
#include "misc.h"

static int tree_error, len;

/****************************************************************************
  Node arena. Nodes are never freed one by one: trees share subtrees, and
  kc builds them for a single run. TreeFreeAll releases all of them.
*****************************************************************************/

#define BlockSize 4096

typedef struct Block {
  struct Block *next;
  int           used;
  TreeNode      nodes[BlockSize];
} Block;

static Block *blocks=NULL;

static TreeNode *NewNode(TreeTag tag) {

  TreeNode *t;

  if ((blocks==NULL) || (blocks->used==BlockSize)) {
    Block *b=malloc(sizeof(Block));
    if (b==NULL) {
      fprintf(stderr, "ERROR: Out of memory.\nAborting kc - sorry.\n");
      exit(-1);
    }
    b->next=blocks;
    b->used=0;
    blocks=b;
  }
  t=&blocks->nodes[blocks->used++];
  memset(t, 0, sizeof(*t));
  t->tag=tag;
  return t;
} /* NewNode */

void TreeFreeAll(void) {

  while (blocks!=NULL) {
    Block *next=blocks->next;
    free(blocks);
    blocks=next;
  }
} /* TreeFreeAll */

/****************************************************************************
  Unreduced node builders.
*****************************************************************************/

static Tree Op(TreeOper op, Tree left, Tree right) {

  TreeNode *t=NewNode(Oper);

  t->data.oper.op=op;
  t->data.oper.left=left;
  t->data.oper.right=right;
  return t;
} /* Op */

static Tree Const(double value) {

  TreeNode *t=NewNode(Konst);

  t->data.value=value;
  return t;
} /* Const */

static Tree Fun(Function fun, Tree arg) {

  TreeNode *t=NewNode(Func);

  t->data.func.fun=fun;
  t->data.func.right=arg;
  return t;
} /* Fun */

static Tree Neg(Tree arg) {

  TreeNode *t=NewNode(NegSign);

  t->data.right=arg;
  return t;
} /* Neg */

static int IsConst(Tree t, double value) {

  return (t->tag==Konst) && (t->data.value==value);
} /* IsConst */

static int SameVar(Tree a, Tree b) {

  return (a->tag==Var) && (b->tag==Var)
    && (strcmp(a->data.name, b->data.name)==0);
} /* SameVar */

/****************************************************************************
  Reduce returns a simpler tree with the same value. It reduces the
  subtrees, then applies the first matching rule below to the node. The
  notation f, g stands for any expression, K for a constant.
*****************************************************************************/

static Tree Reduce(Tree t) {

  Tree l, r;

  if (t==NULL) {
    fprintf(stderr, "ERROR in TreeReduce: NULL pointer.\nAborting kc - sorry.\n");
    exit(-1);
  }
  switch (t->tag) {
  case Oper:
    l=Reduce(t->data.oper.left);
    r=Reduce(t->data.oper.right);
    switch (t->data.oper.op) {
    case Add:
      if ((l->tag==Konst) && (r->tag==Konst))      /* K1+K2 -> K3 */
	return Const(l->data.value+r->data.value);
      if (IsConst(r, 0.0)) return l;               /* f+0 -> f */
      if (IsConst(l, 0.0)) return r;               /* 0+f -> f */
      if (r->tag==NegSign)                         /* f+(-g) -> f-g */
	return Op(Sub, l, r->data.right);
      break;
    case Sub:
      if ((l->tag==Konst) && (r->tag==Konst))      /* K1-K2 -> K3 */
	return Const(l->data.value-r->data.value);
      if (IsConst(r, 0.0)) return l;               /* f-0 -> f */
      if (IsConst(l, 0.0)) return Neg(r);          /* 0-f -> -f */
      if (SameVar(l, r)) return Const(0.0);        /* x-x -> 0 */
      if (r->tag==NegSign)                         /* f-(-g) -> f+g */
	return Op(Add, l, r->data.right);
      if (l->tag==NegSign)                         /* (-f)-g -> -(f+g) */
	return Neg(Op(Add, l->data.right, r));
      break;
    case Mul:
      if ((l->tag==Konst) && (r->tag==Konst))      /* K1*K2 -> K3 */
	return Const(l->data.value*r->data.value);
      if (IsConst(r, 0.0)) return Const(0.0);      /* f*0 -> 0 */
      if (IsConst(r, 1.0)) return l;               /* f*1 -> f */
      if (IsConst(l, 0.0)) return Const(0.0);      /* 0*f -> 0 */
      if (IsConst(l, 1.0)) return r;               /* 1*f -> f */
      break;
    case Div:
      if ((l->tag==Konst) && (r->tag==Konst))      /* K1/K2 -> K3 */
	return Const(l->data.value/r->data.value);
      if (IsConst(r, 1.0)) return l;               /* f/1 -> f */
      if (IsConst(r, 0.0)) return Const(0.0);      /* f/0 -> 0 (sic) */
      if (SameVar(l, r)) return Const(1.0);        /* x/x -> 1 */
      if (IsConst(l, 0.0)) return Const(0.0);      /* 0/f -> 0 */
      break;
    case Pow:
      if ((l->tag==Konst) && (r->tag==Konst))      /* K1^K2 -> K3 */
	return Const(pow(l->data.value, r->data.value));
      if (IsConst(r, 1.0)) return l;               /* f^1 -> f */
      if (IsConst(r, 0.0)) return Const(1.0);      /* f^0 -> 1 */
      break;
    }
    return Op(t->data.oper.op, l, r);
  case NegSign:
    r=Reduce(t->data.right);
    if (r->tag==NegSign) return r->data.right;     /* -(-f) -> f */
    if (r->tag==Konst) return Const(-r->data.value); /* -(K) -> (-K) */
    return Neg(r);
  case Func:
    return Fun(t->data.func.fun, Reduce(t->data.func.right));
  default:
    return t;
  }
} /* Reduce */


static void LineBreak(FILE *out, int mode) {

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

int TreeGetError(void) {

  return tree_error;
} /* TreeGetError */

Tree TreeConst(double value) {

  tree_error=0;
  return Const(value);
} /* TreeConst */

Tree TreeVar(const char *name) {

  TreeNode *t=NewNode(Var);

  (void) strcpy(t->data.name, name);
  tree_error=0;
  return t;
} /* TreeVar */

static Tree Reduced(Tree t) {

  t=Reduce(t);
  tree_error=0;
  return t;
} /* Reduced */

Tree TreeAdd(Tree a, Tree b) { return Reduced(Op(Add, a, b)); }
Tree TreeSub(Tree a, Tree b) { return Reduced(Op(Sub, a, b)); }
Tree TreeMul(Tree a, Tree b) { return Reduced(Op(Mul, a, b)); }
Tree TreeDiv(Tree a, Tree b) { return Reduced(Op(Div, a, b)); }
Tree TreePow(Tree a, Tree b) { return Reduced(Op(Pow, a, b)); }
Tree TreeSign(Tree a) { return Reduced(Neg(a)); }
Tree TreeApplyFunc(Function fun, Tree a) { return Reduced(Fun(fun, a)); }

/* Replaces the variable name by value. Like the original implementation,
   it reduces at every level of the tree. */
Tree TreeSubstTree(Tree t, const char *name, Tree value) {

  Tree res;

  if (t==NULL) {
    tree_error=1;
    return NULL;
  }
  switch (t->tag) {
  case Var:
    res=(strcmp(t->data.name, name)==0) ? value : t;
    break;
  case Oper:
    res=Op(t->data.oper.op, TreeSubstTree(t->data.oper.left, name, value),
	   TreeSubstTree(t->data.oper.right, name, value));
    break;
  case Func:
    res=Fun(t->data.func.fun, TreeSubstTree(t->data.func.right, name, value));
    break;
  case NegSign:
    res=Neg(TreeSubstTree(t->data.right, name, value));
    break;
  default:
    res=t;
    break;
  }
  return Reduced(res);
} /* TreeSubstTree */

Tree TreeSubstVar(Tree t, const char *name, double value) {

  return TreeSubstTree(t, name, Const(value));
} /* TreeSubstVar */

/* Derivative of t with respect to name. Every level is reduced once. */
Tree TreeDerive(Tree t, const char *name) {

  Tree res=NULL, f, g, u, du;

  tree_error=0;
  if (t==NULL)
    return Const(0.0);
  switch (t->tag) {
  case Konst:
    res=Const(0.0);
    break;
  case Var:
    res=Const((strcmp(t->data.name, name)==0) ? 1.0 : 0.0);
    break;
  case Oper:
    f=t->data.oper.left;
    g=t->data.oper.right;
    switch (t->data.oper.op) {
    case Add:
      res=Op(Add, TreeDerive(f, name), TreeDerive(g, name));
      break;
    case Sub:
      res=Op(Sub, TreeDerive(f, name), TreeDerive(g, name));
      break;
    case Mul:                    /* f'*g + f*g' */
      res=Op(Add, Op(Mul, TreeDerive(f, name), g),
	     Op(Mul, f, TreeDerive(g, name)));
      break;
    case Div:                    /* (g*f' - f*g')/g^2 */
      res=Op(Div, Op(Sub, Op(Mul, g, TreeDerive(f, name)),
		     Op(Mul, f, TreeDerive(g, name))), Op(Pow, g, Const(2.0)));
      break;
    case Pow:                    /* g*f^(g-1)*f' + f^g*ln(f)*g' */
      res=TreeAdd(TreeMul(TreeMul(TreePow(f, TreeSub(g, Const(1.0))),
				  TreeDerive(f, name)), g),
		  TreeMul(TreeMul(TreePow(f, g), TreeApplyFunc(Ln, f)),
			  TreeDerive(g, name)));
      break;
    }
    break;
  case Func:
    u=t->data.func.right;
    du=TreeDerive(u, name);
    switch (t->data.func.fun) {
    case Exp:   res=Op(Mul, du, Fun(Exp, u)); break;
    case Sin:   res=Op(Mul, du, Fun(Cos, u)); break;
    case Cos:   res=Op(Mul, du, Neg(Fun(Sin, u))); break;
    case Sinh:  res=Op(Mul, du, Fun(Cosh, u)); break;
    case Cosh:  res=Op(Mul, du, Fun(Sinh, u)); break;
    case Ln:    res=Op(Mul, du, Op(Div, Const(1.0), u)); break;
    case Log:   res=Op(Mul, du, Op(Div, Const(0.4342944819), u)); break;
    case Tan:   res=Op(Div, du, Op(Pow, Fun(Cos, u), Const(2.0))); break;
    case Tanh:  res=Op(Div, du, Op(Pow, Fun(Cosh, u), Const(2.0))); break;
    case Atan:  res=Op(Div, du, Op(Add, Const(1.0), Op(Pow, u, Const(2.0))));
      break;
    case Atanh: res=Op(Div, du, Op(Sub, Const(1.0), Op(Pow, u, Const(2.0))));
      break;
    case Asinh: res=Op(Div, du, Op(Pow, Op(Add, Op(Pow, u, Const(2.0)),
					  Const(1.0)), Const(0.5)));
      break;
    case Acosh: res=Op(Div, du, Op(Pow, Op(Sub, Op(Pow, u, Const(2.0)),
					  Const(1.0)), Const(0.5)));
      break;
    case Asin:  res=Op(Div, du, Op(Pow, Op(Sub, Const(1.0),
					  Op(Pow, u, Const(2.0))), Const(0.5)));
      break;
    case Acos:  res=Op(Div, Neg(du), Op(Pow, Op(Sub, Const(1.0),
					       Op(Pow, u, Const(2.0))), Const(0.5)));
      break;
    }
    break;
  case NegSign:
    res=Neg(TreeDerive(t->data.right, name));
    break;
  }
  return Reduce(res);
} /* TreeDerive */

double TreeEval(Tree t) {
 
  double c1, c2, res=0.0;

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
	  res=asinh(res);
	  break;
	case Acosh:
	  res=acosh(res);
	  break;
	case Atanh:
	  res=atanh(res);
	  break;
	} /* switch */
      } /* if */
      break;
    } /* switch */
  } /* else */
 return res;
} /* TreeEval */

/* The following procedure is 'internal' */

static void TreePrintXtra(Tree t, int mode, FILE *output) {

/* Mode:  1 - (HP) Fortran 77 */
/*	  2 - (HP) Pascal     */
/*        3 - ANSI C          */

  int i;

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
      default:
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
      if (fmod(t->data.value, 1.0)==0.0 && fabs(t->data.value)<=INT_MAX) {
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
