/***************************************************************************
   Misc. functions.

   CopyWrong 1993-1995 by
   Kenneth Geisshirt (kneth@fatou.ruc.dk)
   Department of Life Sciences and Chemistry
   Roskilde University
   P.O. Box 260
   4000 Roskilde
   Denmark
   
   See kc.tex for details.

   Last updated: 15 February 1995 by KG
****************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "config.h"
#include "tableman.h"
#include "misc.h"

void write0(char *str) {

  printf(str);
  printf("\n");
}

void write1(char *str, double r0) {

  printf(str, r0);
  printf("\n");
}

void write2(char *str, double r0, double r1) {
  
  printf(str, r0, r1);
  printf("\n");
}

void write3(char *str, double r0, double r1, double r2) {

  printf(str, r0, r1, r2);
  printf("\n");
}

void GetAndPrintConst(char *name, char *text, int type, double def, 
		      FILE *output, int mode) {

  double temp;
  char   equal[10], line_end[10], prefix[10];

  switch (mode) {
  case 1: /* F77 */
    strcpy(equal, "=");
    strcpy(line_end, "");
    strcpy(prefix, "      ");
    break;
  case 2: /* Pascal */
    strcpy(equal, ":=");
    strcpy(line_end, ";");
    strcpy(prefix, "");
    break;
  case 3: /* C */
    strcpy(equal, "=");
    strcpy(line_end, ";");
    strcpy(prefix, "");
    break;
  } /* switch */
  
  temp=GetConstant(name);
  if (GetError()!=NoError)
    temp=def;
  switch (type) {
  case 0: /* integer */
    fprintf(output, "%s%s%s%d%s\n", prefix, text, equal, (int)temp, line_end);
    break;
  case 1: /* real */
    fprintf(output, "%s%s%s%e%s\n", prefix, text, equal, temp, line_end);
    break;
  } /* switch */
} /* GetAndPrintConst */


int Fact(int n) {

  int i, temp = 1;

  for(i=2; i<=n; i++) 
    temp*=i;
  return temp;
}

int ipower(int x, int n) {

  int i, temp;

  temp=x;
  for(i=2; i<=n; i++) 
    temp*=x;
  return temp;
}

char *StringAlloc(void) {

  char *temp;

  temp=(char *)calloc(STRING_LENGTH, sizeof(char));
  if (temp==NULL)
    fprintf(stderr, "Error in StringAlloc: No space for string.\n");
  return (temp);
} /* StringAlloc */

void StringFree(char *str) {

  free((MALLOCTYPE *)str);
} /* StringFree */
