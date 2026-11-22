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

static int expr_error, len;

/****************************************************************************
  Node arena. Nodes are never freed one by one: trees share subtrees, and
  kc builds them for a single run. expr_free_all releases all of them.
*****************************************************************************/

#define BlockSize 4096

typedef struct Block {
  struct Block *next;
  int           used;
  expr_node_t   nodes[BlockSize];
} Block;

static Block *blocks=NULL;

static expr_node_t *NewNode(expr_tag_t tag) {

  expr_node_t *t;

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

void expr_free_all(void) {

  while (blocks!=NULL) {
    Block *next=blocks->next;
    free(blocks);
    blocks=next;
  }
} /* expr_free_all */

/****************************************************************************
  Unreduced node builders.
*****************************************************************************/

static expr_t Op(expr_op_t op, expr_t left, expr_t right) {

  expr_node_t *t=NewNode(Oper);

  t->data.oper.op=op;
  t->data.oper.left=left;
  t->data.oper.right=right;
  return t;
} /* Op */

/* Constants are never -0, which C prints as -0.000000e+00. */
static expr_t Const(double value) {

  expr_node_t *t=NewNode(Konst);

  t->data.value=(value==0.0) ? 0.0 : value;
  return t;
} /* Const */

static expr_t Fun(expr_func_t fun, expr_t arg) {

  expr_node_t *t=NewNode(Func);

  t->data.func.fun=fun;
  t->data.func.right=arg;
  return t;
} /* Fun */

static expr_t Neg(expr_t arg) {

  expr_node_t *t=NewNode(NegSign);

  t->data.right=arg;
  return t;
} /* Neg */

static int IsConst(expr_t t, double value) {

  return (t->tag==Konst) && (t->data.value==value);
} /* IsConst */

static int SameVar(expr_t a, expr_t b) {

  return (a->tag==Var) && (b->tag==Var)
    && (strcmp(a->data.name, b->data.name)==0);
} /* SameVar */

/****************************************************************************
  Reduction rules. Each function takes reduced operands and returns a
  reduced tree, applying the first matching rule. Every tree outside this
  file is reduced, since nodes are only built here and every exported
  routine goes through these rules; so a new node needs only its own
  rules, and reduced subtrees are never walked again. The notation f, g
  stands for any expression, K for a constant.
*****************************************************************************/

static expr_t MulR(expr_t l, expr_t r);
static expr_t SubR(expr_t l, expr_t r);

/* Is t of the form K*f? MulR keeps constant factors on the left. */
static int IsScaled(expr_t t) {

  return (t->tag==Oper) && (t->data.oper.op==Mul)
    && (t->data.oper.left->tag==Konst);
} /* IsScaled */

static double Factor(expr_t t) {

  return t->data.oper.left->data.value;
} /* Factor */

/* Is t of the form -f, K or K*f with K<0? */
static int IsNegative(expr_t t) {

  return (t->tag==NegSign) || ((t->tag==Konst) && (t->data.value<0.0))
    || (IsScaled(t) && (Factor(t)<0.0));
} /* IsNegative */

static expr_t NegR(expr_t r) {

  if (r->tag==NegSign) return r->data.right;     /* -(-f) -> f */
  if (r->tag==Konst) return Const(-r->data.value); /* -(K) -> (-K) */
  if (IsScaled(r))                               /* -(K*f) -> (-K)*f */
    return MulR(Const(-Factor(r)), r->data.oper.right);
  return Neg(r);
} /* NegR */

static expr_t AddR(expr_t l, expr_t r) {

  if ((l->tag==Konst) && (r->tag==Konst))        /* K1+K2 -> K3 */
    return Const(l->data.value+r->data.value);
  if (IsConst(r, 0.0)) return l;                 /* f+0 -> f */
  if (IsConst(l, 0.0)) return r;                 /* 0+f -> f */
  if (IsNegative(r))                   /* f+(-g) -> f-g, see IsNegative */
    return SubR(l, NegR(r));
  return Op(Add, l, r);
} /* AddR */

static expr_t SubR(expr_t l, expr_t r) {

  if ((l->tag==Konst) && (r->tag==Konst))        /* K1-K2 -> K3 */
    return Const(l->data.value-r->data.value);
  if (IsConst(r, 0.0)) return l;                 /* f-0 -> f */
  if (IsConst(l, 0.0)) return NegR(r);           /* 0-f -> -f */
  if (SameVar(l, r)) return Const(0.0);          /* x-x -> 0 */
  if (IsNegative(r))                   /* f-(-g) -> f+g, see IsNegative */
    return AddR(l, NegR(r));
  if (l->tag==NegSign)                           /* (-f)-g -> -(f+g) */
    return NegR(AddR(l->data.right, r));
  return Op(Sub, l, r);
} /* SubR */

static expr_t MulR(expr_t l, expr_t r) {

  if ((l->tag==Konst) && (r->tag==Konst))        /* K1*K2 -> K3 */
    return Const(l->data.value*r->data.value);
  if (IsConst(l, 0.0) || IsConst(r, 0.0))        /* 0*f, f*0 -> 0 */
    return Const(0.0);
  if (r->tag==Konst) return MulR(r, l);          /* f*K -> K*f */
  if (IsConst(l, 1.0)) return r;                 /* 1*f -> f */
  if (IsConst(l, -1.0)) return NegR(r);          /* -1*f -> -f */
  if ((l->tag==Konst) && IsScaled(r))            /* K1*(K2*f) -> K3*f */
    return MulR(Const(l->data.value*Factor(r)), r->data.oper.right);
  if (IsScaled(l))                               /* (K*f)*g -> K*(f*g) */
    return MulR(l->data.oper.left, MulR(l->data.oper.right, r));
  if (IsScaled(r))                               /* f*(K*g) -> K*(f*g) */
    return MulR(r->data.oper.left, MulR(l, r->data.oper.right));
  if (l->tag==NegSign)                           /* (-f)*g -> -(f*g) */
    return NegR(MulR(l->data.right, r));
  if (r->tag==NegSign)                           /* f*(-g) -> -(f*g) */
    return NegR(MulR(l, r->data.right));
  return Op(Mul, l, r);
} /* MulR */

static expr_t DivR(expr_t l, expr_t r) {

  if ((l->tag==Konst) && (r->tag==Konst))        /* K1/K2 -> K3 */
    return Const(l->data.value/r->data.value);
  if (IsConst(r, 1.0)) return l;                 /* f/1 -> f */
  if (IsConst(r, 0.0)) return Const(0.0);        /* f/0 -> 0 (sic) */
  if (SameVar(l, r)) return Const(1.0);          /* x/x -> 1 */
  if (IsConst(l, 0.0)) return Const(0.0);        /* 0/f -> 0 */
  if (l->tag==NegSign)                           /* (-f)/g -> -(f/g) */
    return NegR(DivR(l->data.right, r));
  if (IsScaled(l))                               /* (K*f)/g -> K*(f/g) */
    return MulR(l->data.oper.left, DivR(l->data.oper.right, r));
  return Op(Div, l, r);
} /* DivR */

static expr_t PowR(expr_t l, expr_t r) {

  if ((l->tag==Konst) && (r->tag==Konst))        /* K1^K2 -> K3 */
    return Const(pow(l->data.value, r->data.value));
  if (IsConst(r, 1.0)) return l;                 /* f^1 -> f */
  if (IsConst(r, 0.0)) return Const(1.0);        /* f^0 -> 1 */
  return Op(Pow, l, r);
} /* PowR */

static expr_t OpR(expr_op_t op, expr_t l, expr_t r) {

  switch (op) {
  case Add: return AddR(l, r);
  case Sub: return SubR(l, r);
  case Mul: return MulR(l, r);
  case Div: return DivR(l, r);
  case Pow: return PowR(l, r);
  }
  return Op(op, l, r);
} /* OpR */


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

int expr_get_error(void) {

  return expr_error;
} /* expr_get_error */

expr_t expr_const(double value) {

  return Const(value);
} /* expr_const */

expr_t expr_var(const char *name) {

  expr_node_t *t=NewNode(Var);

  (void) strcpy(t->data.name, name);
  return t;
} /* expr_var */

expr_t expr_add(expr_t a, expr_t b) { return AddR(a, b); }
expr_t expr_sub(expr_t a, expr_t b) { return SubR(a, b); }
expr_t expr_mul(expr_t a, expr_t b) { return MulR(a, b); }
expr_t expr_div(expr_t a, expr_t b) { return DivR(a, b); }
expr_t expr_pow(expr_t a, expr_t b) { return PowR(a, b); }
expr_t expr_neg(expr_t a) { return NegR(a); }
expr_t expr_apply_func(expr_func_t fun, expr_t a) { return Fun(fun, a); }

/* Replaces the variable name by value. Subtrees without name are
   shared, not copied. */
expr_t expr_subst_expr(expr_t t, const char *name, expr_t value) {

  expr_t l, r;

  if (t==NULL)
    return NULL;
  switch (t->tag) {
  case Var:
    return (strcmp(t->data.name, name)==0) ? value : t;
  case Oper:
    l=expr_subst_expr(t->data.oper.left, name, value);
    r=expr_subst_expr(t->data.oper.right, name, value);
    if ((l==t->data.oper.left) && (r==t->data.oper.right)) return t;
    return OpR(t->data.oper.op, l, r);
  case Func:
    r=expr_subst_expr(t->data.func.right, name, value);
    return (r==t->data.func.right) ? t : Fun(t->data.func.fun, r);
  case NegSign:
    r=expr_subst_expr(t->data.right, name, value);
    return (r==t->data.right) ? t : NegR(r);
  default:
    return t;
  }
} /* expr_subst_expr */

expr_t expr_subst_var(expr_t t, const char *name, double value) {

  return expr_subst_expr(t, name, Const(value));
} /* expr_subst_var */

static expr_t Sq(expr_t t) { return PowR(t, Const(2.0)); }
static expr_t Sqrt(expr_t t) { return PowR(t, Const(0.5)); }

/* Derivative of t with respect to name. */
expr_t expr_derive(expr_t t, const char *name) {

  expr_t f, g, u, du;

  if (t==NULL)
    return Const(0.0);
  switch (t->tag) {
  case Var:
    return Const((strcmp(t->data.name, name)==0) ? 1.0 : 0.0);
  case Oper:
    f=t->data.oper.left;
    g=t->data.oper.right;
    switch (t->data.oper.op) {
    case Add:
      return AddR(expr_derive(f, name), expr_derive(g, name));
    case Sub:
      return SubR(expr_derive(f, name), expr_derive(g, name));
    case Mul:                    /* f'*g + f*g' */
      return AddR(MulR(expr_derive(f, name), g), MulR(f, expr_derive(g, name)));
    case Div:                    /* (g*f' - f*g')/g^2 */
      return DivR(SubR(MulR(g, expr_derive(f, name)),
		       MulR(f, expr_derive(g, name))), Sq(g));
    case Pow:                    /* g*f^(g-1)*f' + f^g*ln(f)*g' */
      return AddR(MulR(MulR(PowR(f, SubR(g, Const(1.0))),
			    expr_derive(f, name)), g),
		  MulR(MulR(PowR(f, g), Fun(Ln, f)), expr_derive(g, name)));
    }
    break;
  case Func:
    u=t->data.func.right;
    du=expr_derive(u, name);
    switch (t->data.func.fun) {
    case Exp:   return MulR(du, t);
    case Sin:   return MulR(du, Fun(Cos, u));
    case Cos:   return MulR(du, NegR(Fun(Sin, u)));
    case Sinh:  return MulR(du, Fun(Cosh, u));
    case Cosh:  return MulR(du, Fun(Sinh, u));
    case Ln:    return MulR(du, DivR(Const(1.0), u));
    case Log:   return MulR(du, DivR(Const(0.4342944819), u));
    case Tan:   return DivR(du, Sq(Fun(Cos, u)));
    case Tanh:  return DivR(du, Sq(Fun(Cosh, u)));
    case Atan:  return DivR(du, AddR(Const(1.0), Sq(u)));
    case Atanh: return DivR(du, SubR(Const(1.0), Sq(u)));
    case Asinh: return DivR(du, Sqrt(AddR(Sq(u), Const(1.0))));
    case Acosh: return DivR(du, Sqrt(SubR(Sq(u), Const(1.0))));
    case Asin:  return DivR(du, Sqrt(SubR(Const(1.0), Sq(u))));
    case Acos:  return DivR(NegR(du), Sqrt(SubR(Const(1.0), Sq(u))));
    }
    break;
  case NegSign:
    return NegR(expr_derive(t->data.right, name));
  default:
    break;
  }
  return Const(0.0);
} /* expr_derive */

double expr_eval(expr_t t) {
 
  double c1, c2, res=0.0;

  if (t==NULL) {
    expr_error=EXPR_NO_EXPR;
    return 0.0;
  }
  switch (t->tag) {
  case Var:
    expr_error=EXPR_NO_EVAL;
    res=0.0;
    break;
  case Konst:
    res=t->data.value;
    expr_error=EXPR_OK;
    break;
  case Oper: 
    c1=expr_eval(t->data.oper.left);
    if (expr_error==EXPR_OK) {
      c2=expr_eval(t->data.oper.right);
      if (expr_error==EXPR_OK) {
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
      } /* if */
    } /* if */
    break;
  case NegSign:
    res=-expr_eval(t->data.right);
    break;
  case Func:
    res=expr_eval(t->data.func.right);
    if (expr_error==EXPR_OK) {
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
 return res;
} /* expr_eval */

/* The following procedure is 'internal' */

static void expr_print_xtra(expr_t t, int mode, FILE *output) {

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
      expr_print_xtra(t->data.oper.left, mode, output);
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
      expr_print_xtra(t->data.oper.right, mode, output);
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
      expr_print_xtra(t->data.func.right, mode, output);
      fprintf(output, ")");
      len++;
      LineBreak(output, mode);
      break; /* Func */
    case NegSign:
      fprintf(output, "(-(");
      len+=2;
      LineBreak(output, mode);
      expr_print_xtra(t->data.right, mode, output);
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
	expr_print_xtra(t->data.oper.left, mode, output);
	fprintf(output, "+");
	len++;
	LineBreak(output, mode);
	expr_print_xtra(t->data.oper.right, mode, output);
	break;
      case Sub: 
	expr_print_xtra(t->data.oper.left, mode, output);
	fprintf(output, "-"); 
	len++;
	LineBreak(output, mode);
	expr_print_xtra(t->data.oper.right, mode, output);
	break;
      case Div: 
	expr_print_xtra(t->data.oper.left, mode, output);
	fprintf(output, "/"); 
	len++;
	LineBreak(output, mode);
	expr_print_xtra(t->data.oper.right, mode, output);
	break;
      case Mul: 
	expr_print_xtra(t->data.oper.left, mode, output);
	fprintf(output, "*"); 
	len++;
	LineBreak(output, mode);
	expr_print_xtra(t->data.oper.right, mode, output);
	break;
      case Pow: 
	if (t->data.oper.right->tag==Konst) {
	  if (ceil(t->data.oper.right->data.value)==floor(t->data.oper.right->data.value)) {
	    expr_print_xtra(t->data.oper.left, mode, output);
	    for(i=1; i<((int)(t->data.oper.right->data.value)); i++) {
	      fprintf(output, "*");
	      expr_print_xtra(t->data.oper.left, mode, output);
	      len++;
	      LineBreak(output, mode);
	    }
	  } else {
	    fprintf(output, "EXP("); /* Pascal has no power operator!! */ 
	    len+=4;
	    LineBreak(output, mode);
	    expr_print_xtra(t->data.oper.right, mode, output);
	    fprintf(output, "*LOG(");
	    len+=5;
	    LineBreak(output, mode);
	    expr_print_xtra(t->data.oper.left, mode, output);
	    fprintf(output, "))");
	    len+=2;
	    LineBreak(output, mode);
	  } /* if ... else ... */
	} else {
	  fprintf(output, "EXP("); /* Pascal has no power operator!! */ 
	  len+=4;
	  LineBreak(output, mode);
	  expr_print_xtra(t->data.oper.right, mode, output);
	  fprintf(output, "*LOG(");
	  len+=5;
	  LineBreak(output, mode);
	  expr_print_xtra(t->data.oper.left, mode, output);
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
      expr_print_xtra(t->data.func.right, mode, output);
      fprintf(output, ")");
      len++;
      LineBreak(output, mode);
      break; /* Func */
    case NegSign:
      fprintf(output, "(-(");
      len+=2;
      LineBreak(output, mode);
      expr_print_xtra(t->data.right, mode, output);
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
	expr_print_xtra(t->data.oper.left, mode, output);
	fprintf(output, "/");
	len++;
	LineBreak(output, mode);
	expr_print_xtra(t->data.oper.right, mode, output);
	break;
      case Mul: 
	expr_print_xtra(t->data.oper.left, mode, output);
	fprintf(output, "*");
	len++;
	LineBreak(output, mode);
	expr_print_xtra(t->data.oper.right, mode, output);
	break;
      case Add: 
	expr_print_xtra(t->data.oper.left, mode, output);
	fprintf(output, "+");
	len++;
	LineBreak(output, mode);
	expr_print_xtra(t->data.oper.right, mode, output);
	break;
      case Sub:
	expr_print_xtra(t->data.oper.left, mode, output);
	fprintf(output, "-");
	len++;
	LineBreak(output, mode);
	expr_print_xtra(t->data.oper.right, mode, output);
	break;
      case Pow:
	if (t->data.oper.right->tag==Konst) {
	  if (ceil(t->data.oper.right->data.value)==floor(t->data.oper.right->data.value)) {
	    expr_print_xtra(t->data.oper.left, mode, output);
	    for(i=1; i<((int)(ceil(t->data.oper.right->data.value))); i++) {
	      fprintf(output, "*");
	      len++;
	      LineBreak(output, mode);
	      expr_print_xtra(t->data.oper.left, mode, output);
	    } /* for i */
	  } else {
	    fprintf(output, "pow(");
	    len+=4;
	    LineBreak(output, mode);
	    expr_print_xtra(t->data.oper.left, mode, output);
	    fprintf(output, ",");
	    len++;
	    LineBreak(output, mode);
	    expr_print_xtra(t->data.oper.right, mode, output);
	    fprintf(output, ")");
	    len++;
	    LineBreak(output, mode);
	  } /* if ... else ... */
	} else {
	  fprintf(output, "pow(");
	  len+=4;
	  LineBreak(output, mode);
	  expr_print_xtra(t->data.oper.left, mode, output);
	  fprintf(output, ",");
	  len++;
	  LineBreak(output, mode);
	  expr_print_xtra(t->data.oper.right, mode, output);
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
      expr_print_xtra(t->data.func.right, mode, output);
      fprintf(output, ")");
      len++;
      LineBreak(output, mode);
      break; /* Func */
    case NegSign:
      fprintf(output, "(-(");
      len+=2;
      LineBreak(output, mode);
      expr_print_xtra(t->data.right, mode, output);
      fprintf(output, "))");
      len++;
      LineBreak(output, mode);
      break; /* NegSign */
    }; /* switch tag */
    break; /* mode 3 */   
  }; /* switch */
} /* expr_print_xtra */


void expr_print(expr_t t, int mode, FILE *output) {

  len=10;
  expr_print_xtra(t, mode, output);
} /* expr_print */
