/************************************************************************** 
   Misc. routines for kc.

   CopyWrong 1993-1995 by
   Kenneth Geisshirt (kneth@fatou.ruc.dk)
   Department of Life Sciences and Chemistry
   Roskilde University
   4000 Roskilde
   Denmark

   See kc.tex for details.

   Last updated: 15 February 1995 by KG
****************************************************************************/

#ifndef _MISC_
#define _MISC_

#include <stdio.h>

int         verbose;         /* 1==verbose on */
extern void write0(char *);
extern void write1(char *, double);
extern void write2(char *, double, double);
extern void write3(char *, double, double, double);
extern void GetAndPrintConst(char *, char *, int, double, FILE *, int);
extern int  Fact(int);
extern int  ipower(int, int);
extern char *StringAlloc(void);
extern void StringFree(char *);
#endif
