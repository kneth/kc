/**************************************************************************
  This code generator is an experimental one. It came from a discuss
  with Dimitrios Pananakis, Dundee University, School of Biomedical
  Engineering. In general there are no error checks, and activation 
  energy, temperature power and reference temperature are not used
  and generated.

  CopyWrong 1994 by
  Kenneth Geisshirt (kneth@osc.kiku.dk)
  Department of Theoretical Chemistry
  H.C. Orsted Institute
  Universitetsparken 5
  2100 Copenhagen
  Denmark

  Last updated: 11 April 1994
****************************************************************************/

#include <stdio.h>
#include <malloc.h>

#include "config.h"
#include "symbmath.h"
#include "tableman.h"

void Pananakis(FILE *code) {

  char *name, *rename;
  int i, j, l, react_no;
  double temp, charge, coeff;
  Tree temp_tree;
  
  name=(char *)malloc(sizeof(char));
  rename=(char *)malloc(sizeof(char));
  
  for(i=1; i<=NoOfSpec(); i++) {
    GetSpecNo(i, name, &charge);
    fprintf(code, "REACTION(SET=%d,TERMS=%d)\n", i, NoSpecInReacs(name, charge));
    for(j=1; j<=NoOfReact(); j++) {
      react_no=GetReactNo(j-1);
      GetSpecNo(i, name, &charge);
      if (IsSpecInReact(react_no, name, charge, &coeff)==1) {
	switch (GetReactKind(react_no)) {
	case uni:
	  temp_tree=TreeCreate();
	  GetRateConst(react_no, uni, 1, temp_tree);
	  temp=TreeEval(temp_tree);
	  if (TreeGetError()==NoEval)
	    fprintf(stderr, "Unable to compute rate constant in reaction %d.\n", react_no);
	  TreeKill(temp_tree);
	  fprintf(code, "%e", -coeff*temp);
	  for(l=1; l<=NoOfSpec(); l++) {
	    GetSpecNo(l, name, &charge);
	    if (IsSpecInReact(react_no, name, charge, &coeff)==1) {
	      if (coeff>0.0)
		fprintf(code, ",%e", coeff);
	      else
		fprintf(code, ",0");
	    } else
	      fprintf(code, ",0");
	  } /* for l */
	  break;
	case bi:
	  temp_tree=TreeCreate();
	  GetRateConst(react_no, bi, 1, temp_tree);
	  temp=TreeEval(temp_tree);
	  if (TreeGetError()==NoEval)
	    fprintf(stderr, "Unable to compute rate constant (>) in reaction %d.\n", react_no);
	  TreeKill(temp_tree);
	  fprintf(code, "%e", -coeff*temp);
	  for(l=1; l<=NoOfSpec(); l++) {
	    GetSpecNo(l, name, &charge);
	    if (IsSpecInReact(react_no, name, charge, &coeff)==1) {
	      if (coeff>0.0)
		fprintf(code, ",%e", coeff);
	      else
		fprintf(code, ",0");
	    } else
	      fprintf(code, ",0");
	  } /* for l */
	  fprintf(code, "\n");
	  temp_tree=TreeCreate();
	  GetRateConst(react_no, bi, 2, temp_tree);
	  temp=TreeEval(temp_tree);
	  if (TreeGetError()==NoEval)
	    fprintf(stderr, "Unable to compute rate constant (<) in reaction %d.\n", react_no);
	  TreeKill(temp_tree);
	  fprintf(code, "%e", -coeff*temp);
	  for(l=1; l<=NoOfSpec(); l++) {
	    GetSpecNo(l, name, &charge);
	    if (IsSpecInReact(react_no, name, charge, &coeff)==1) {
	      if (coeff<0.0)
		fprintf(code, ",%e", -coeff);
	      else
		fprintf(code, ",0");
	    } else
	      fprintf(code, ",0");
	  } /* for l */
	  break;
	} /* switch */
	fprintf(code, "\n");
      }
    } /* for j */
  } /* for i */
} /* Pananakis */
