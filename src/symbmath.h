/****************************************************************************
  Symbolic mathematics in ANSI-C.

  CopyWrong 1992-1994 by
  Kenneth Geisshirt (kneth@osc.kiku.dk)
  Department of Theoretical Chemistry
  H.C. Orsted Institute
  Universitetsparken 5
  2100 Copenhagen
  Denmark

  See kc.tex for details.

  Expression trees are immutable: every operation returns a new tree and
  shares the subtrees of its arguments. All nodes are allocated from one
  arena, which TreeFreeAll releases; there is no per-tree free.
****************************************************************************/

#ifndef _SYMBMATH_
#define _SYMBMATH_
#include <stdio.h>
#include "config.h"

/* Status of the last TreeEval, returned by TreeGetError. */
#define TreeNoError 1
#define NoEval  2
#define NoTree  3

typedef enum {Konst, Var, Oper, Func, NegSign} expr_tag_t;
typedef enum {Add, Sub, Mul, Div, Pow} expr_op_t;
typedef enum {Exp, Sin, Cos, Tan, Ln, Log, Cosh, Sinh, Tanh, Asin,
		Acos, Atan, Acosh, Asinh, Atanh} expr_func_t;

struct expr_node {
	  expr_tag_t tag;
	  union {
	    double value;    /* constant */
	    char name[STRING_LENGTH];       /* variable */
            struct {
	      expr_op_t op;
	      const struct expr_node *left, *right;
	    } oper;          /* operation */
	    struct {
	      expr_func_t fun;
	      const struct expr_node *right;
	    } func;
	    const struct expr_node *right;      /* negative sign */
          } data;
};
typedef struct expr_node expr_node_t;
typedef const expr_node_t *expr_t;

extern int    TreeGetError(void);
extern expr_t TreeConst(double);
extern expr_t TreeVar(const char *);
extern expr_t TreeAdd(expr_t, expr_t);
extern expr_t TreeSub(expr_t, expr_t);
extern expr_t TreeMul(expr_t, expr_t);
extern expr_t TreeDiv(expr_t, expr_t);
extern expr_t TreePow(expr_t, expr_t);
extern expr_t TreeSign(expr_t);
extern expr_t TreeApplyFunc(expr_func_t, expr_t);
extern expr_t TreeSubstVar(expr_t, const char *, double);
extern expr_t TreeSubstTree(expr_t, const char *, expr_t);
extern expr_t TreeDerive(expr_t, const char *);
extern double TreeEval(expr_t);
extern void   TreePrint(expr_t, int, FILE *);
extern void   TreeFreeAll(void);
#endif
