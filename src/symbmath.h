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

#define TreeNoError 1
#define NoEval  2
#define NoTree  3

typedef enum {Konst, Var, Oper, Func, NegSign} TreeTag;
typedef enum {Add, Sub, Mul, Div, Pow} TreeOper;
typedef enum {Exp, Sin, Cos, Tan, Ln, Log, Cosh, Sinh, Tanh, Asin,
		Acos, Atan, Acosh, Asinh, Atanh} Function;

struct TreeCell {
	  TreeTag   tag;
	  union {
	    double value;    /* constant */
	    char name[STRING_LENGTH];       /* variable */
            struct {
	      TreeOper op;
	      const struct TreeCell *left, *right;
	    } oper;          /* operation */
	    struct {
	      Function fun;
	      const struct TreeCell *right;
	    } func;
	    const struct TreeCell *right;      /* negative sign */
          } data;
};
typedef struct TreeCell TreeNode;
typedef const TreeNode *Tree;

extern int    TreeGetError(void);
extern Tree   TreeConst(double);
extern Tree   TreeVar(const char *);
extern Tree   TreeAdd(Tree, Tree);
extern Tree   TreeSub(Tree, Tree);
extern Tree   TreeMul(Tree, Tree);
extern Tree   TreeDiv(Tree, Tree);
extern Tree   TreePow(Tree, Tree);
extern Tree   TreeSign(Tree);
extern Tree   TreeApplyFunc(Function, Tree);
extern Tree   TreeSubstVar(Tree, const char *, double);
extern Tree   TreeSubstTree(Tree, const char *, Tree);
extern Tree   TreeDerive(Tree, const char *);
extern double TreeEval(Tree);
extern void   TreePrint(Tree, int, FILE *);
extern void   TreeFreeAll(void);
#endif
