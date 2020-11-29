/* DSCML - Dalimil Snita's Chemical Meta Langauge
   CopyWrong by Kenneth Geisshirt, 1992, 1993

   See kc.tex for details
*/

#include <stdlib.h>
#include <malloc.h>
#include <stdio.h>
#include "config.h"
#include "symbmath.h"
#include "tableman.h"

void PrintCharge(FILE *code, double charge) {

/* This routine prints a proper string for the charge */

  int i;
  
  if (charge!=0.0) {
    fprintf(code, "(");
    if (charge==MAXFLOAT)
      fprintf(code, "rad");
    else if (charge>0.0) 
      for(i=1; i<=(int)charge; i++)
        fprintf(code, "+");
    else if (charge<0.0)
      for(i=1; i<=-(int)charge; i++)
         fprintf(code, "-");
    fprintf(code, ")");
  }; /* if */
} /* PrintCharge */

void DSCML(FILE *code) {

  int    i, j, finished, no;
  char   *name;
  double charge, temp, coeff;
  Tree   tmp;
  
  name=(char *)malloc(sizeof(char));
  fprintf(code, "%e\t\t\t- number of components\n\n", (double)NoOfSpec());
  for(i=1; i<=NoOfSpec(); i++) {
    GetSpecNo(i, name, &charge);
    fprintf(code, "%s\t\t\t- name of %d. component\n", name, i);
    fprintf(code, "%e\t\t\t- charge number\n", charge);
    temp=GetSpecConst(name, charge, "M");
    if (GetError()==NotFound)
      fprintf(stderr, "Molar mass is not defined for %s(%d).\n", name, (int)charge);
    else
      fprintf(code, "%e\t\t\t- molar mass\n", temp);
    temp=GetSpecConst(name, charge, "D");
    if (GetError()==NotFound)
      fprintf(stderr, "Diffusivity is not defined for %s(%d).\n", name, (int)charge);
    else
      fprintf(code, "%e\t\t\t- diffusivity\n", temp);
    temp=GetBeginConc(name, charge);
    fprintf(code, "%e\t\t\t- initial concentration\n\n", temp);
  }; /* for i */
  fprintf(code, "%d\t\t\t- number of reactions\n\n", NoOfReact());
  j=1; 
  for(i=1; i<=NoOfReact(); i++) {
    no=GetReactNo(i-1);
    switch (GetReactKind(no)) {
    case uni: 
      finished=GetFirstSpecA(no, name, &charge, &coeff, 0);
      fprintf(code, "'&d. reaction\t\t\t- name of the %d. reaction\n", j, j);
      tmp=TreeCreate();
      GetRateConst(no, uni, 1, tmp);
      TreePrint(tmp, 1, code);
      TreeKill(tmp);
      fprintf(code, "\t\t\t- kinetic constant\n");
      while (finished==0) {
	fprintf(code, "%e\t\t\t- %s(", coeff, name);
	PrintCharge(code, charge);
	fprintf(code, ")\n");
	finished=GetNextSpecA(name, &charge, &coeff, 0);
      }; /* while */
      fprintf(code, "\n");
      j++;
      break; /* uni */
    case bi:  
      finished=GetFirstSpecA(no, name, &charge, &coeff, 0); 
      fprintf(code, "'&d. reaction\t\t\t- name of the %d. reaction\n", j, j);
      tmp=TreeCreate();
      GetRateConst(no, bi, 1, tmp);
      TreePrint(tmp, 2, code);
      TreeKill(tmp);
      fprintf(code, "\t\t\t- kinetic constant\n", temp);
      while (finished==0) {
	fprintf(code, "%e\t\t\t- %s(", coeff, name);
	PrintCharge(code, charge);
	fprintf(code, ")\n");
	finished=GetNextSpecA(name, &charge, &coeff, 0);
      }; /* while */
      fprintf(code, "\n");
      j++;
      fprintf(code, "'&d. reaction\t\t\t- name of the %d. reaction\n", j, j);
      tmp=TreeCreate();
      GetRateConst(no, bi, 2, tmp);
      TreePrint(tmp, 2, code);
      TreeKill(tmp);
      fprintf(code, "\t\t\t- kinetic constant\n", temp);
      finished=GetFirstSpecA(no, name, &charge, &coeff, 1);
      while (finished==0) {
	fprintf(code, "%e\t\t\t- %s(", coeff, name);
	PrintCharge(code, charge);
	fprintf(code, ")\n");
	finished=GetNextSpecA(name, &charge, &coeff, 1);
      }; /* while */
      fprintf(code, "\n");
      j++;
      break; /* bi */
    case equi:  
      fprintf(stderr, "Equilibriums are not supported.\n");
      break; /* equi */
    }; /* switch */
  }; /* for i */
} /* DSCML */

