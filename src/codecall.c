/*************************************************************************
  CodeGenCall - callee of all code generators.

  CopyWrong 1993-1995 by 
  Kenneth Geisshirt (kneth@fatou.ruc.dk)
  Department of Life Sciences and Chemistry
  Roskilde University
  P.O. Box 260
  4000 Roskilde
  Denmark

  See kc.tex for details.

  Last updated: 9 May 1995 by KN
**************************************************************************/


#include <stdio.h>
#include "config.h"
#include "tableman.h"

/* import all code generators */
#include "waves.h"
#include "kncont.h"
#include "kgode.h"
#include "conis.h"

void CodeGenCall(int mode) {

/************************************************************************** 
  modes:
  1 - Mixed - Different small programs (idea by K. Nielsen)
  2 - Waves - per1d/KGadi (K. Geisshirt)
  3 - KGode/kci (K. Geisshirt and K. Nielsen)
  4 - Output to a continuation program by I. Schreiber (K. Nielsen)
  5 - KnCont (K. Nielsen)
***************************************************************************/

  FILE *code, *code_h, *code_c, *code_ini, *code_ini2;
  int  intgr;

  switch (mode) {
  case 1:
    Mixed();
    break; 
  case 2:
    code_h=fopen("model.h", "w");
    code_c=fopen("model.c", "w");
    code_ini=fopen("in.dat", "w");
    Waves(code_h, code_c, code_ini);
    fclose(code_ini);
    fclose(code_c);
    fclose(code_h);
    break; 
  case 3:
    code_c=fopen("model.c", "w");
    code_h=fopen("model.h", "w");
    KGode(code_c, code_h);
    fclose(code_c);
    fclose(code_h);
    break;
  case 4:
    code_c=fopen("model.f", "w");
    ConIS(code_c);
    fclose(code_c);
    break;
  case 5:
    code_c=fopen("kcm5proc.p", "w");
    code_h=fopen("kcm5const.p", "w");
    KnCont(code_c, code_h);
    fclose(code_c);
    fclose(code_h);
    break;
  default:
    fprintf(stderr, "CodeGenCall: mode %d is not supported. Run kc -h for help.\n", mode);
  }; /* switch */
} /* CodeGenCall */
