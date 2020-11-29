/*************************************************************************** 
  Main program for kc. It is actually just parsing arguments, gives control 
  to the parser, and then calls the code generator..

  (C) Copyright 1992-1996 by 
  Kenneth Geisshirt (kneth@fatou.ruc.dk)           Keld Neilsen (kn@kiku.dk)
  Dept. of Life Sciences and Chemistry        Dept. of Theoretical Chemistry
  Roskilde University                               University of Copenhagen
  P.O. Box 260                                          Universitetsparken 5
  4000 Roskilde                                         2100 Copenhagen East
  Denmark                                                            Denmark
  
  See kc.tex for details.

  Last updated: 4 January 1996 by KG
***************************************************************************/


#define DefaultMode 1

#include "config.h"
#include <stdio.h>
#ifdef _USE_GARBAGE_COL_
#  include <gc.h>
#else
#  if _MALLOC_DEBUG_
#     include "malloc.h"
#  else
#    include <stdlib.h>
#  endif
#endif
#include "codecall.h"
#include "misc.h"

int    mode;

#include "parser.c"

int main(int argc, char *argv[]) {

  int    j, c;
  extern char *optarg;

#ifdef _MALLOC_DEBUG_
  union dbmalloptarg m;
  dbmallopt(MALLOC_ERRFILE, m);
#endif

#ifdef _PLATFORM_HPUX_
  mallopt(M_MXFAST, sizeof(Tree));  
#endif

  mode=DefaultMode;
  printf("kc v%s, (C) Copyright by Kenneth Geisshirt and Keld Nielsen, 1992-1996\n", VERSION);
  verbose=1;                                             /* printing  */
  while ((c=getopt(argc, argv, "vhdm:"))!=EOF)   
    switch (c) {
    case 'h': 
      (void) printf("\nOptions:\n");
      (void) printf("  -h : this text.\n");
      (void) printf("  -v : verbose off\n");
      (void) printf("  -mx: mode\n");
      (void) printf("       1 = Calculation of various properties\n");
      (void) printf("       2 = KGadi\n");
      (void) printf("       3 = kci\n");
      (void) printf("       4 = Continuation program by I.Schreiber\n");
      (void) printf("       5 = KnCont\n");
      (void) printf("E-mail: kneth@fatou.ruc.dk or kn@kiku.dk\n\n");
      (void) printf("This program comes with absolutely no warrenty!\n");
      exit(-1);
      break;
    case 'm': 
      j=0;
      mode=0;
      while ((optarg[j]!=' ') && (optarg[j]!='\0')) {
	mode=mode*10+(int)(optarg[j]-'0');
	j++;
      };
      break;
    case 'v':
      verbose=0;
      break;
    default: 
      (void) fprintf(stderr, "Warning: Option %1s ignored.\n", argv[i][1]);
      break;
    } /* switch */
  (void) yyparse();
  CodeGenCall(mode);
} 
