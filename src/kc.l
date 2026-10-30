%option noyywrap
%option nounput
%{
static void CommentC(void);
static void CommentPas(void);
static void num(void);
%}
Digit       [0-9]
Letter      [A-Za-z]
Sign        "-"|"+"
alphabet    ({Digit}|{Letter})
%%
[ \t] ;
\n   lineno++;
"/*" { CommentC(); }
"(*" { CommentPas(); }

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
{Letter}({Digit}|{Letter})* { (void) strcpy(yylval.name, yytext); return names; } 
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
static void CommentC(void)
{ char c;
  int finished;

  finished=0;
  while (finished==0) {
    while (input()!='*') /* nothing */;
    if ((c=input())=='/') finished=1;
  } /* while */
} /* Comment */       

static void CommentPas(void)
{ char c;
  int finished;

  finished=0;
  while (finished==0) {
    while (input()!='*') /* nothing */;
    if ((c=input())==')') finished=1;
  } /* while */
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
