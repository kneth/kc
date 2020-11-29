/* AutoCont - a code generator for kc and CONT using 
   power law and mass-action law.
   CopyWrong by Kenneth Geisshirt, 1992

   See kc.tex for details.
   Version 0.00

   Only for experimental use!!!
*/

#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <time.h>

static int line_len;

void LineBreaker(FILE *out) {

  if (line_len>=60) { 
    fprintf(out, "\n     &"); 
    line_len=0; 
  };
} /* LineBreaker */

void AutoCont(FILE *code) {

  int i, j;                               /* counters                 */
  int  in_use[ReactSize];                 /* any right to left reactions? */
  char *name, *rename;                    /* formal and actual names  */
  double charge, coeff, value;            /* misc. values             */
  int react_no;                           /* reaction number (source) */
  int finished;                           /* boolean stop variable    */
  int line_len;                           /* #chars written on line   */
  time_t timer;

  name=malloc(sizeof(char));
  rename=malloc(sizeof(char));
  timer=time(&timer);	    
	    
/* The actual code generator */

	    fprintf(code, "      SUBROUTINE MODEL(NDIM,NVAR,N,T,X,F,G)\n");
	    fprintf(code, "C *********************************************\n");
	    fprintf(code, "C %s", ctime(&timer));
	    fprintf(code, "C *********************************************\n");
	    fprintf(code, "      IMPLICIT REAL *8(A-H,O-Z)\n");
	    for(i=0; i<=NoOfSpec()-1; i++) {
	      GetSpecNo(i+1, name, &charge);
	      RenameSpec(rename, name, charge);
	      fprintf(code, "      REAL %s\n", rename);
            }; /* for */
	    fprintf(code, "      DIMENSION X(NDIM),F(NDIM),G(NDIM,NVAR)\n");
	    fprintf(code, "      COMMON/FIXP/PAR(20)\n");
	    fprintf(code, "      COMMON/VARP/ALPHA,BETA,ARG,PER\n");
	    fprintf(code, "\nC Make aliases\n");
	    for(i=0; i<=NoOfSpec()-1; i++) {
	      GetSpecNo(i+1, name, &charge);
	      RenameSpec(rename, name, charge);
	      fprintf(code, "      %s = X(%d)\n", rename, i);
            }; /* for */
	    for(i=1; i<=NoOfReact(); i++) {
              fprintf(code, "C Rate of reaction no. ");  
	      react_no=GetReactNo(i-1);  /* not pretty, but needed */
	      fprintf(code, "%d\n", react_no);
	      switch (GetReactKind(react_no)) {
		case uni  : line_len=20;
			    fprintf(code, "      F(%d) = ", i-1);
	                    fprintf(code, "%e", GetRateConst(react_no, uni, 1)); 
			    finished=GetFirstSpecA(react_no, name, &charge, &coeff);
                            while (finished==1) {
		              if (coeff>0.0) {
				RenameSpec(rename, name, charge);
			        LineBreaker(code);	
				fprintf(code, "*(%s**%e)", rename, coeff);
				line_len=line_len+strlen(rename)+13;
                              }; /* if */
		              finished=GetNextSpecA(name, &charge, &coeff);
                            }; /* while */
			    fprintf(code, "\n");
			    break;
		case bi   : /* Left to right reaction */
			    line_len=20;
			    fprintf(code, "      F(%d) = ", i-1);
			    fprintf(code, "%e", GetRateConst(react_no, bi, 1)); 
	                    finished=GetFirstSpecA(react_no, name, &charge, &coeff);
                            while (finished==1) {
		              if (coeff>0.0) {
			        RenameSpec(rename, name, charge);
				LineBreaker(code);
				fprintf(code, "*(%s**%e)", rename, coeff);
				line_len=line_len+strlen(rename)+13;
                              }; /* if */
		              finished=GetNextSpecA(name, &charge, &coeff);
                            }; /* while */
                            /* Right to left reaction */
			    fprintf(code, "-%e", GetRateConst(react_no, bi, 2));
			    line_len=line_len+10;
	                    finished=GetFirstSpecA(react_no, name, &charge, &coeff);
                            while (finished==1) {
		              if (coeff<0.0) {
				RenameSpec(rename, name, charge);
				LineBreaker(code);
				fprintf(code, "*(%s**%e)", rename, -coeff);
                                line_len=line_len+strlen(rename)+13;
			      }; /* if */
		              finished=GetNextSpecA(name, &charge, &coeff);
                            }; /* while */
			    fprintf(code, "\n");
			    break;
		case equi : /* not implemented */
			    fprintf(stderr, "Warning: Equilibriums are not implemented as feature - yet\n"); 
			    break;
	      }; /* switch */
            }; /* for */
            fprintf(code, "C *******************************************\n");
	    fprintf(code, "      IF (N.EQ.NDIM) RETURN\n");
	    fprintf(code, "C *******************************************\n");
	    fprintf(code, "C Jacobian:\n");
	    fprintf(code, "      RETURN\n");
	    fprintf(code, "      END\n");
} /* AutoCont */
