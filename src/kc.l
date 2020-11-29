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
"K"      return K;    
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
CommentC() 
{ char c;
  int finished;

  finished=0;
  while (finished==0) {
    while (input()!='*') /* nothing */;
    if ((c=input())=='/') finished=1;
  } /* while */
} /* Comment */       

CommentPas() 
{ char c;
  int finished;

  finished=0;
  while (finished==0) {
    while (input()!='*') /* nothing */;
    if ((c=input())==')') finished=1;
  } /* while */
} /* Comment */       


incl() {

  int i, start;

  i=8; 
  while (yytext[i]==' ') i++;
  start=i;
  while (yytext[i]!='\0') { yytext[i-start]=yytext[i]; i++; }
  yytext[i]='\0';
} /* incl */

num() {
  
  float temp;
  int i;

  for(i=0; i<strlen(yytext); i++)
    if ((yytext[i]=='L') || (yytext[i]=='l'))
      yytext[i]='e';
#ifdef _PLATFORM_AIX_
  sscanf((char *)yytext, "%e", &temp);
#else
  sscanf(yytext, "%e", &temp);
#endif
  value=(double)temp;
} /* num */

int yyerror(char *s) {

  fprintf(stderr, "Line %d: %s at input text %s\n", lineno, s, yytext);
} /* yyerror */
