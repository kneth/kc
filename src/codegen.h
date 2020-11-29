/**************************************************************************** 
  CodeGen - code generator for kc. The routines in this library should be
  called by the actual code generator.

  CopyWrong 1992-1995 by
  Kenneth Geisshirt (kneth@fatou.ruc.dk)
  Department of Life Sciences and Chemistry
  Roskilde University
  P.O. Box 260
  4000 Roskilde
  Denmark

  See kc.tex for details.

  Last updated: 11 May 1995 by KN
******************************************************************************/

#ifndef _CODEGEN_
#define _CODEGEN_

#include "symbmath.h"
#include "tableman.h"

/* Global variables which can be used by the code generators. */
Tree *v;
Tree *con;
Tree **jacobi;
Tree ***hess;
Tree ****keld;
Tree *rfw;
Tree *rrv;
Tree **rfwds;
Tree **rrvds;


/* exported routines */
extern void InitCodeGenVar(int, int, int);
extern void GenerateRateExpr(void);
extern void GenerateJacobi(void);
extern void GenerateHessian(void);
extern void EvalJacoby(double *, double **);
extern void GenerateKeldian(void);
extern void GenerateDiffReac(void);
#endif
