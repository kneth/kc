/*************************************************************************** 
  Stoc - code generator to kc.
  The code generator is viewing reactions as a stocastic
  process, see [1].

  CopyWrong, 1993, 1994 by
  Kenneth Geisshirt (kneth@osc.kiku.dk)
  Department of Theoretical Chemistry
  H.C. Orsted Institute
  Universitetsparken 5
  2100 Copenhagen
  Denmark

  See kc.tex for details.

  References:
  [1] A general method for numerically simulating the stochastic time 
      evolution of coupled chemical reactions.
      D. T. Gillespie, J. Comp. Phys., vol 22, pp. 403-434 (1976).

  Last updated: 4 April 1994
*****************************************************************************/

#include "config.h"
#include "tableman.h"
#include "symbmath.h"
#include "misc.h"

#include <stdio.h>
#include <malloc.h>
#include <time.h>
#include <math.h>

void StocCode(FILE *code) {

  char      *name, *rename;
  time_t    timer;
  double    volume;
  int       i, j, react_no, finished, dyn;
  Tree      temp_tree;
  double    coeff, charge, temp;
  int       fac, sfac;
  double    rate_const[2*MaxReact];

  name=(char *)malloc(sizeof(char));
  rename=(char *)malloc(sizeof(char));

  timer=time(&timer);
  fprintf(code, "/******************************************************\n");
  fprintf(code, "  Warning: This file was generated by kc v%s\n", VERSION);
  fprintf(code, "  CopyWrong by Kenneth Geisshirt (kneth@osc.kiku.dk)\n");
  fprintf(code, "  %s*/\n", ctime(&timer));
  fprintf(code, "*******************************************************/\n");

  volume=GetConstant("vol");
  if (GetError()==NotFound)
    fprintf(stderr, "STOC: Volume (vol) is not speficied.\n");
  
  dyn=0;
  for(i=1; i<=NoOfReact(); i++) {
    react_no=GetReactNo(i-1);
    switch (GetReactKind(react_no)) {
    case uni:
      if (GetRateKind(react_no, uni, 1)==1) {
	dyn++;
	temp_tree=TreeCreate();
	GetRateConst(react_no, uni, 1, temp_tree);
	temp=TreeEval(temp_tree);
	if (TreeGetError()==NoEval)
	  fprintf(stderr, "STOC: Unable to evaluate rate constant for reaction %d\n", react_no);
	TreeKill(temp_tree);
	rate_const[2*(i-1)]=pow(volume, SumCoeff(1))*temp;
	finished=GetFirstSpecA(react_no, name, &charge, &coeff, 0);
	while (finished==1) {
	  constraint=IsSpecInConstraint(name, charge);
	  if (constraint==0)
	    rate_const[2*(i-1)]/=(double)Fact((int)coeff);
	  finished=GetNextSpecA(name, &charge, &coeff, 0);
	} /* while */
      } /* if */
      break;
    case bi:
      if (GetRateKind(react_no, bi, 1)==1) {
	dyn++;
	temp_tree=TreeCreate();
	GetRateConst(react_no, bi, 1, temp_tree);
	temp=TreeEval(temp_tree);
	if (TreeGetError()==NoEval)
	  fprintf(stderr, "STOC: Unable to evaluate rate constant for reaction %d\n", react_no);
	TreeKill(temp_tree);
	rate_const[2*(i-1)]=pow(volume, SumCoeff(1))*temp;
	finished=GetFirstSpecA(react_no, name, &charge, &coeff, 0);
	while (finished==1) {
	  constraint=IsSpecInConstraint(name, charge);
	  if (constraint==0)
	    rate_const[2*(i-1)]/=(double)Fact((int)coeff);
	  finished=GetNextSpecA(name, &charge, &coeff, 0);
	} /* while */
      } /* if */
      if (GetRateKind(react_no, bi, 2)==1) {
	dyn++;
	temp_tree=TreeCreate();
	GetRateConst(react_no, bi, 2, temp_tree);
	temp=TreeEval(temp_tree);
	if (TreeGetError()==NoEval)
	  fprintf(stderr, "STOC: Unable to evaluate rate constant for reaction %d\n", react_no);
	rate_const[2*(i-1)+1]=pow(volume, SumCoeff(2))*temp;
	finished=GetFirstSpecA(react_no, name, &charge, &coeff, 1);
	while (finished==1) {
	  constraint=IsSpecInConstraint(name, charge);
	  if (constraint==0)
	    rate_const[2*(i-1)+1].value/=(double)Fact((int)coeff);
	  finished=GetNextSpecA(name, &charge, &coeff, 1);
	} /* while */
      } /* if */
      break;
    case equi:
      fprintf(stderr, "STOC: Cannot handle equilibriums.\n");
      break;
    } /* switch */
  } /* for i */

  fprintf(code, "double A[%d];\n", dyn);

  fprintf(code, "void CalcA(void) {\n");
  dyn=-1;
  for(i=1; i<=NoOfReact(); i++) {
    react_no=GetReactNo(i-1);
    switch (GetReactKind(react_no)) {
    case uni:
      if (GetRateKind(react_no, uni, 1)==1) {
	dyn++;
	fprintf(code, "A[%d]=%e", dyn, rate_const[2*(i-1)]);
	finished=GetFirstSpecA(react_no, name, &charge, &coeff, 0);
	while (finished==1) {
	  RenameSpec(rename, name, charge);
	  constraint=IsSpecInConstraint(name, charge);
	  if (constraint==0) {
	    for(j=1; j<=(int)coeff; j++)
	      fprintf(code, "*(%s-%d)", rename, j-1);
	    fprintf(code, "/%e", (double)Fact((int)coeff));
	  }
	  finished=GetNextSpecA(name, &charge, &coeff, 0);
	} /* while */
      } /* if */
      break;
    case bi:
      if (GetRateKind(react_no, bi, 1)==1) {
	dyn++;
	fprintf(code, "A[%d]=%e", dyn, rate_const[2*(i-1)]);
	finished=GetFirstSpecA(react_no, name, &charge, &coeff, 0);
	while (finished==1) {
	  RenameSpec(rename, name, charge);
	  constraint=IsSpecInConstraint(name, charge);
	  if (constraint==0) {
	    for(j=1; j<=(int)coeff; j++)
	      fprintf(code, "*(%s-%d)", rename, j-1);
	    fprintf(code, "/%e", (double)Fact((int)coeff));
	  }
	  finished=GetNextSpecA(name, &charge, &coeff, 0);
	} /* while */
      } /* if */
     if (GetRateKind(react_no, bi, 2)==1) {
	dyn++;
	fprintf(code, "A[%d]=%e", dyn, rate_const[2*(i-1)+1]);
	finished=GetFirstSpecA(react_no, name, &charge, &coeff, 1);
	while (finished==1) {
	  RenameSpec(rename, name, charge);
	  constraint=IsSpecInConstraint(name, charge);
	  if (constraint==0) {
	    for(j=1; j<=(int)coeff; j++)
	      fprintf(code, "*(%s-%d)", rename, j-1);
	    fprintf(code, "/%e", (double)Fact((int)coeff));
	  }
	  finished=GetNextSpecA(name, &charge, &coeff, 1);
	} /* while */
      } /* if */
      break;
    case equi:
      fprintf(stderr, "STOC: Cannot handle equilibriums.\n");
      break;
    } /* switch */
  } /* for i */
  
    

