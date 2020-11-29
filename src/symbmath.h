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

  Last updated: 27 September 1994
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
	      struct TreeCell *left, *right;
	    } oper;          /* operation */
	    struct {
	      Function fun;
	      struct TreeCell *right;
	    } func;
	    struct TreeCell *right;      /* negative sign */
          } data;
};
typedef struct TreeCell TreeNode;
typedef TreeNode *Tree;

int           tree_error, len, left, right;

extern int    TreeGetError();
extern Tree   TreeCreate(void);
extern void   TreeAdd(Tree, Tree);
extern void   TreeSub(Tree, Tree);
extern void   TreeMul(Tree, Tree);
extern void   TreeDiv(Tree, Tree);
extern void   TreePow(Tree, Tree);
extern void   TreeSign(Tree);
extern void   TreeAssignConst(Tree, double);
extern void   TreeAssignVar(Tree, char *);
extern void   TreeSubstVar(Tree, char *, double);
extern void   TreeDerive(Tree, Tree, char *);
extern double TreeEval(Tree);
extern void   TreePrint(Tree, int, FILE *);
extern void   TreeCpy(Tree, Tree);
extern void   TreeKill(Tree);
extern void   TreeSubstTree(Tree, char *, Tree);
extern void   TreeApplyFunc(Tree *, Function);  
#endif
