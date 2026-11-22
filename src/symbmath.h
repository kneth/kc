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
  arena, which expr_free_all releases; there is no per-tree free.
****************************************************************************/

#ifndef _SYMBMATH_
#define _SYMBMATH_
#include <stdio.h>
#include "config.h"

/* Status of the last expr_eval, returned by expr_get_error. */
#define EXPR_OK      1
#define EXPR_NO_EVAL 2
#define EXPR_NO_EXPR 3

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

extern int    expr_get_error(void);
extern expr_t expr_const(double);
extern expr_t expr_var(const char *);
extern expr_t expr_add(expr_t, expr_t);
extern expr_t expr_sub(expr_t, expr_t);
extern expr_t expr_mul(expr_t, expr_t);
extern expr_t expr_div(expr_t, expr_t);
extern expr_t expr_pow(expr_t, expr_t);
extern expr_t expr_neg(expr_t);
extern expr_t expr_apply_func(expr_func_t, expr_t);
extern expr_t expr_subst_var(expr_t, const char *, double);
extern expr_t expr_subst_expr(expr_t, const char *, expr_t);
extern expr_t expr_derive(expr_t, const char *);
extern double expr_eval(expr_t);
extern void   expr_print(expr_t, int, FILE *);
extern void   expr_free_all(void);
#endif
