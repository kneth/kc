%option noyywrap
%option nounput
%{
static void Comment(int close);
static void num(void);
%}
Digit       [0-9]
Letter      [A-Za-z]
Sign        "-"|"+"
alphabet    ({Digit}|{Letter})
%%
[ \t] ;
\n   lineno++;
"/*" { Comment('/'); }
"(*" { Comment(')'); }

"exp"    return fun_exp;
"log"    return fun_log;
"ln"     return fun_ln;
"sin"    return fun_sin;
"cos"    return fun_cos;
"tan"    return fun_tan;
"sinh"   return fun_sinh;
"cosh"   return fun_cosh;
"tanh"   return fun_tanh;
"asin"   return fun_asin;
"acos"   return fun_acos;
"atan"   return fun_atan;
"asinh"  return fun_asinh;
"acosh"  return fun_acosh;
"atanh"  return fun_atanh;

(({Digit}+)|({Digit}*\.{Digit}+)([eELl][-+]?{Digit}+)?) { num(); yylval.dval=value; return numbers; }
{Letter}({Digit}|{Letter})* { if (yyleng>=STRING_LENGTH) {
                                fprintf(stderr, "Line %d: name %s is longer than %d characters\n",
                                        lineno, yytext, STRING_LENGTH-1);
                                exit(EXIT_FAILURE);
                              } /* if */
                              (void) strcpy(yylval.name, yytext); return names; }
"#parameter" return param;
"#print" return print;
"'"      return prime;
"("      return leftpar;
")"      return rightpar;
"(0)"    return time0;
"->"     return oneway;
"<->"    return twoways;
"<=>"    return twoways;
"="      return equal;
"."      return radical;
"k>"     { strcpy(yylval.name, "k"); return leftarr; }
"k<"     { strcpy(yylval.name, "k"); return rightarr; }
"v>"     { strcpy(yylval.name, "v"); return leftarr; }
"v<"     { strcpy(yylval.name, "v"); return rightarr; }
"+"      return plus;
"-"      return minus;
"**"     return powop;
"^"      return powop;
"*"      return multi;
"/"      return pdiv;
"["      return leftconc;
"]"      return rightconc;
","      return comma;
";"      return semicolon;
":"      return colon;
"\""     return quotation;

%%
/* Skips a comment up to and including the closing '*' followed by CLOSE. */
static void Comment(int close)
{ int c, prev=0, start=lineno;

  while ((c=input())!=EOF && c!=0) {
    if (c=='\n')
      lineno++;
    if ((prev=='*') && (c==close))
      return;
    prev=c;
  } /* while */
  fprintf(stderr, "Line %d: comment is not closed\n", start);
  exit(EXIT_FAILURE);
} /* Comment */

static void num(void) {
  
  float temp;
  size_t i;

  for(i=0; i<strlen(yytext); i++)
    if ((yytext[i]=='L') || (yytext[i]=='l'))
      yytext[i]='e';
  sscanf(yytext, "%e", &temp);
  value=(double)temp;
} /* num */

int yyerror(char *s) {

  fprintf(stderr, "Line %d: %s at input text %s\n", lineno, s, yytext);
  return 0;
} /* yyerror */
