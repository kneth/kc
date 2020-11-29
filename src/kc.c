
# line 1 "kc.y"
 
/********************************************************* 
  Parser to kc.
  CopyWrong 1992-1995 by 
  Kenneth Geisshirt (kneth@fatou.ruc.dk)
  Department of Life Sciences and Chemistry
  Roskilde University
  P.O. Box 260
  4000 Roskilde
  Denmark


  See kc.tex for details.

  Last updated: 17 February 1994
**********************************************************/

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include "symbmath.h"
#include "tableman.h"
#include "codecall.h"

typedef enum {no, ini, ord} Conc; 

typedef struct compound {
  char   name[STRING_LENGTH];
  double charge;
  Conc   concs;
} compound;

typedef struct Strnum {
  int    flag;
  char   name[STRING_LENGTH];
  double numb;
} Strnum;

double     coeff, charge, value, temp, temp1, temp2;
char       flag, name[STRING_LENGTH], string[STRING_LENGTH];
int        i, j, side, lineno=1;
Tree       tmp;

# line 46 "kc.y"
typedef union  {
  double      dval;
  char        oper;
  char        name[STRING_LENGTH];
  compound    comp;
  char        flag;
  Tree        tree;
  Function    func;
  Strnum      strnum;
} YYSTYPE;
# define names 257
# define leftarr 258
# define rightarr 259
# define numbers 260
# define param 261
# define print 262
# define powop 263
# define semicolon 264
# define quotation 265
# define equal 266
# define colon 267
# define R 268
# define E 269
# define powc 270
# define leftpar 271
# define rightpar 272
# define oneway 273
# define twoways 274
# define plus 275
# define minus 276
# define multi 277
# define pdiv 278
# define comma 279
# define leftconc 280
# define rightconc 281
# define time0 282
# define K 283
# define radical 284
# define V 285
# define prime 286
# define fun_exp 287
# define fun_log 288
# define fun_ln 289
# define fun_sin 290
# define fun_cos 291
# define fun_tan 292
# define fun_sinh 293
# define fun_cosh 294
# define fun_tanh 295
# define fun_asin 296
# define fun_acos 297
# define fun_atan 298
# define fun_asinh 299
# define fun_acosh 300
# define fun_atanh 301
# define UMINUS 302
#define yyclearin yychar = -1
#define yyerrok yyerrflag = 0
extern int yychar;
extern int yyerrflag;
#ifndef YYMAXDEPTH
#define YYMAXDEPTH 150
#endif
YYSTYPE yylval, yyval;
typedef int yytabelem;
#include <stdio.h>
# define YYERRCODE 256
yytabelem yyexca[] ={
	-1, 1,
	0, -1,
	-2, 0,
	};
# define YYNPROD 83
# define YYLAST 375
yytabelem yyact[]={

    45,    31,    24,    46,    98,    99,    31,    33,    42,    71,
    37,   159,    33,    96,    44,   148,    25,    33,    88,    43,
   147,    73,    25,    18,    16,    38,    21,    10,    30,    27,
    49,    51,    50,    52,    53,    54,    55,    56,    57,    58,
    59,    60,    61,    62,    63,    45,    83,    18,    46,    18,
   120,    18,    18,    92,    39,   171,   163,    83,   104,    44,
    81,    82,   146,   134,    43,   102,   103,   101,    18,    79,
    80,    81,    82,   165,   101,    49,    51,    50,    52,    53,
    54,    55,    56,    57,    58,    59,    60,    61,    62,    63,
    83,   117,    70,   162,    69,    83,    67,    66,   145,    64,
    83,   111,    79,    80,    81,    82,   164,    79,    80,    81,
    82,   157,    79,    80,    81,    82,   156,    83,    77,   124,
    11,    65,    83,    10,    12,     9,   121,    83,    83,    79,
    80,    81,    82,   144,    79,    80,    81,    82,   143,    79,
    80,    81,    82,   132,    83,    36,   167,   142,   141,    83,
    84,   140,   139,    75,    83,    76,    79,    80,    81,    82,
   131,    79,    80,    81,    82,   116,    79,    80,    81,    82,
   115,    83,   173,     6,   166,    83,   133,   158,    14,   138,
   137,   122,    83,    79,    80,    81,    82,    79,    80,    81,
    82,   125,    83,   119,    79,    80,    81,    82,    83,   118,
    41,    23,    83,   114,    79,    80,    81,    82,    26,    83,
    79,    80,    81,    82,    79,    80,    81,    82,   112,    83,
    78,    79,    80,    81,    82,    83,    35,    20,    15,    32,
    19,    79,    80,    81,    82,     8,    47,    79,    80,    81,
    82,    17,    29,     7,    85,    86,    22,    34,     5,    28,
    17,    13,     4,     3,     2,     1,    40,   123,    68,    48,
   100,    97,    95,    72,     0,    89,    74,    90,    91,     0,
    93,    94,     0,     0,     0,    22,   129,    87,     0,     0,
   106,   107,   108,   109,   110,     0,     0,     0,     0,   113,
     0,     0,     0,     0,     0,     0,     0,   152,   151,     0,
     0,     0,     0,   105,     0,     0,     0,     0,     0,     0,
     0,     0,     0,     0,     0,     0,   126,   127,   128,     0,
     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     0,   135,   136,     0,     0,     0,     0,     0,     0,
     0,     0,     0,     0,   149,   150,     0,   153,     0,   155,
   130,     0,     0,     0,     0,     0,     0,   160,   161,     0,
     0,     0,     0,     0,   168,   169,   170,     0,     0,     0,
     0,     0,   172,     0,   154 };
yytabelem yypact[]={

 -1000, -1000, -1000,  -137,  -233, -1000, -1000, -1000, -1000,  -231,
 -1000,  -264,  -228,  -229, -1000, -1000,  -270,  -275,  -112,  -254,
 -1000, -1000, -1000,  -213,  -257,  -167,  -143,  -169,  -170, -1000,
  -265,  -112,  -172, -1000,  -174,  -272,  -250, -1000,  -231,  -142,
 -1000,   -44,  -107,  -212,  -212, -1000, -1000,  -275,  -253, -1000,
 -1000, -1000, -1000, -1000, -1000, -1000, -1000, -1000, -1000, -1000,
 -1000, -1000, -1000, -1000,  -212, -1000,  -212,  -212,  -219,  -212,
  -212, -1000, -1000,  -271, -1000,  -208,  -112, -1000, -1000,  -212,
  -212,  -212,  -212,  -212,  -164, -1000,   -54, -1000,  -212,   -61,
  -109,  -114,  -175,   -65,   -71,  -222, -1000,  -134, -1000, -1000,
 -1000,  -142, -1000, -1000, -1000, -1000,  -217,  -217,  -135,  -135,
 -1000,  -145, -1000,   -81, -1000,  -212,  -212,  -212, -1000, -1000,
 -1000, -1000,  -142,  -112, -1000, -1000,  -119,  -136,   -88,  -201,
 -1000,  -212,  -212, -1000,  -110,  -141,  -146, -1000,  -166, -1000,
 -1000,  -204,  -251,  -212,  -212,  -110,  -212,  -112,  -212,  -163,
  -168, -1000, -1000,   -38,  -261,   -38,  -212,  -212,  -171,  -210,
  -173,  -206,  -113,  -212,  -212,  -212, -1000,  -211,   -38,   -38,
   -38,  -212,   -92, -1000 };
yytabelem yypgo[]={

     0,   200,   263,   262,   155,   261,   229,   260,   226,   236,
   153,   259,   256,   255,   254,   253,   252,   251,   248,   243,
   235,   230,   227,   208,   173,   201,   181,   180,   179,   152,
   151,   177,   174,   228 };
yytabelem yyr1[]={

     0,    14,    13,    15,    15,    18,    18,    18,    19,    12,
    12,    21,    21,    22,    22,    20,    23,    23,    16,    16,
    25,    26,    24,    24,    27,    27,    28,    28,    30,    31,
    29,    29,    32,    32,     7,     7,     7,    10,    10,     8,
     2,     2,     3,     3,     3,     5,     5,     4,     4,    17,
    17,    33,    33,    33,     1,     1,     1,     1,     1,     1,
     1,     1,     1,     1,     1,    11,    11,    11,    11,    11,
    11,    11,    11,    11,    11,    11,    11,    11,    11,    11,
     9,     6,     6 };
yytabelem yyr2[]={

     0,     1,     9,     4,     0,     2,     2,     6,     7,     5,
     9,     6,     2,     3,     3,     6,    27,    27,     2,     4,
     1,     1,    18,    11,     6,     2,     6,     2,    13,     1,
    12,     7,     1,     9,     3,     3,     3,     5,     9,     5,
     7,     1,     3,     5,     3,     3,     3,     3,     1,     2,
     4,    11,    15,    11,     5,     7,     7,     7,     7,     7,
     7,     3,     3,     5,     9,     3,     3,     3,     3,     3,
     3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
     7,     1,     3 };
yytabelem yychk[]={

 -1000,   -13,   -14,   -15,   -16,   -18,   -24,   -19,   -20,   262,
   260,   257,   261,   -17,   -24,   -33,   257,    -9,   280,   -21,
   -22,   257,    -9,   -25,   266,   286,   -23,   257,    -9,   -33,
   257,   271,    -6,   282,    -6,    -8,   257,   264,   279,   267,
   -12,    -1,   265,   276,   271,   257,   260,    -9,   -11,   287,
   289,   288,   290,   291,   292,   293,   294,   295,   296,   297,
   298,   299,   300,   301,   266,   264,   266,   266,    -8,   266,
   266,   281,    -2,   271,   -22,   -10,    -4,   260,   264,   275,
   276,   277,   278,   263,   257,    -1,    -1,    -6,   271,    -1,
    -1,    -1,   272,    -1,    -1,    -3,   284,    -5,   275,   276,
    -7,   275,   273,   274,   266,    -8,    -1,    -1,    -1,    -1,
    -1,   265,   272,    -1,   264,   279,   279,   266,   264,   264,
   272,   260,   -26,    -4,   264,   272,    -1,    -1,    -1,   -10,
    -8,   279,   279,   264,   264,    -1,    -1,   -27,   -28,   -29,
   -30,   258,   257,   279,   279,   264,   266,   271,   266,    -1,
    -1,   -29,   -30,    -1,    -8,    -1,   279,   279,   -31,   272,
    -1,    -1,   264,   266,   279,   279,   -32,   259,    -1,    -1,
    -1,   266,    -1,   264 };
yytabelem yydef[]={

     1,    -2,     4,     0,     0,     3,    18,     5,     6,     0,
    20,     0,     0,     2,    19,    49,    81,    81,     0,     0,
    12,    13,    14,     0,     0,     0,     0,     0,     0,    50,
    81,     0,     0,    82,     0,     0,    41,     7,     0,    48,
     8,     0,     0,     0,     0,    61,    62,    81,     0,    65,
    66,    67,    68,    69,    70,    71,    72,    73,    74,    75,
    76,    77,    78,    79,     0,    15,     0,     0,     0,     0,
     0,    80,    39,     0,    11,     0,     0,    47,     9,     0,
     0,     0,     0,     0,     0,    54,     0,    63,     0,     0,
     0,     0,     0,     0,     0,     0,    42,    44,    45,    46,
    21,    48,    34,    35,    36,    37,    55,    56,    57,    58,
    59,     0,    60,     0,    23,     0,     0,     0,    53,    51,
    40,    43,    48,     0,    10,    64,     0,     0,     0,     0,
    38,     0,     0,    52,     0,     0,     0,    22,     0,    25,
    27,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     0,    24,    26,    29,     0,    31,     0,     0,     0,     0,
     0,     0,    32,     0,     0,     0,    30,     0,    28,    16,
    17,     0,     0,    33 };
typedef struct { char *t_name; int t_val; } yytoktype;
#ifndef YYDEBUG
#	define YYDEBUG	0	/* don't allow debugging */
#endif

#if YYDEBUG

char * yyreds[] =
{
	"-no such reduction-",
      "system : /* empty */",
      "system : declars reactions constants",
      "declars : declars declar",
      "declars : /* empty */",
      "declar : constant",
      "declar : parameter",
      "declar : print printlist semicolon",
      "constant : names equal strnumb",
      "strnumb : expr semicolon",
      "strnumb : quotation names quotation semicolon",
      "printlist : printlist comma prn_entry",
      "printlist : prn_entry",
      "prn_entry : names",
      "prn_entry : conc",
      "parameter : param pentry semicolon",
      "pentry : names equal expr comma expr comma expr comma expr comma expr comma expr",
      "pentry : conc equal expr comma expr comma expr comma expr comma expr comma expr",
      "reactions : reaction",
      "reactions : reactions reaction",
      "reaction : numbers",
      "reaction : numbers colon substs kind",
      "reaction : numbers colon substs kind substs semicolon reackonst",
      "reaction : names prime equal expr semicolon",
      "reackonst : powconsts semicolon rateconsts",
      "reackonst : rateconsts",
      "powconsts : powconsts semicolon powconst",
      "powconsts : powconst",
      "powconst : names leftpar subst rightpar equal expr",
      "rateconsts : leftarr equal expr",
      "rateconsts : leftarr equal expr semicolon ropt",
      "rateconsts : names equal expr",
      "ropt : /* empty */",
      "ropt : rightarr equal expr semicolon",
      "kind : oneway",
      "kind : twoways",
      "kind : equal",
      "substs : coeff subst",
      "substs : substs plus coeff subst",
      "subst : names charge",
      "charge : leftpar size rightpar",
      "charge : /* empty */",
      "size : radical",
      "size : sign numbers",
      "size : sign",
      "sign : plus",
      "sign : minus",
      "coeff : numbers",
      "coeff : /* empty */",
      "constants : const",
      "constants : constants const",
      "const : conc timeopt equal expr semicolon",
      "const : names leftpar subst rightpar equal expr semicolon",
      "const : names timeopt equal expr semicolon",
      "expr : minus expr",
      "expr : expr plus expr",
      "expr : expr minus expr",
      "expr : expr multi expr",
      "expr : expr pdiv expr",
      "expr : expr powop expr",
      "expr : leftpar expr rightpar",
      "expr : names",
      "expr : numbers",
      "expr : conc timeopt",
      "expr : function leftpar expr rightpar",
      "function : fun_exp",
      "function : fun_ln",
      "function : fun_log",
      "function : fun_sin",
      "function : fun_cos",
      "function : fun_tan",
      "function : fun_sinh",
      "function : fun_cosh",
      "function : fun_tanh",
      "function : fun_asin",
      "function : fun_acos",
      "function : fun_atan",
      "function : fun_asinh",
      "function : fun_acosh",
      "function : fun_atanh",
      "conc : leftconc subst rightconc",
      "timeopt : /* empty */",
      "timeopt : time0",
};
yytoktype yytoks[] =
{
	"names",	257,
	"leftarr",	258,
	"rightarr",	259,
	"numbers",	260,
	"param",	261,
	"print",	262,
	"powop",	263,
	"semicolon",	264,
	"quotation",	265,
	"equal",	266,
	"colon",	267,
	"R",	268,
	"E",	269,
	"powc",	270,
	"leftpar",	271,
	"rightpar",	272,
	"oneway",	273,
	"twoways",	274,
	"plus",	275,
	"minus",	276,
	"multi",	277,
	"pdiv",	278,
	"comma",	279,
	"leftconc",	280,
	"rightconc",	281,
	"time0",	282,
	"K",	283,
	"radical",	284,
	"V",	285,
	"prime",	286,
	"fun_exp",	287,
	"fun_log",	288,
	"fun_ln",	289,
	"fun_sin",	290,
	"fun_cos",	291,
	"fun_tan",	292,
	"fun_sinh",	293,
	"fun_cosh",	294,
	"fun_tanh",	295,
	"fun_asin",	296,
	"fun_acos",	297,
	"fun_atan",	298,
	"fun_asinh",	299,
	"fun_acosh",	300,
	"fun_atanh",	301,
	"UMINUS",	302,
	"-unknown-",	-1	/* ends search */
};
#endif /* YYDEBUG */

/* @(#)27       1.7.1.3  src/bos/usr/ccs/bin/yacc/yaccpar, cmdlang, bos411, 9432B411a 8/10/94 14:01:53 */
/*
 * COMPONENT_NAME: (CMDLANG) Language Utilities
 *
 * FUNCTIONS: yyparse
 * ORIGINS: 3
 */
/*
** Skeleton parser driver for yacc output
*/

/*
** yacc user known macros and defines
*/
#ifdef YYSPLIT
#   define YYERROR      return(-2)
#else
#   define YYERROR      goto yyerrlab
#endif
#ifdef YACC_MSG
#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE
#endif
#include <nl_types.h>
nl_catd yyusercatd;
#endif
#define YYACCEPT        return(0)
#define YYABORT         return(1)
#ifndef YACC_MSG
#define YYBACKUP( newtoken, newvalue )\
{\
        if ( yychar >= 0 || ( yyr2[ yytmp ] >> 1 ) != 1 )\
        {\
                yyerror( "syntax error - cannot backup" );\
                YYERROR;\
        }\
        yychar = newtoken;\
        yystate = *yyps;\
        yylval = newvalue;\
        goto yynewstate;\
}
#else
#define YYBACKUP( newtoken, newvalue )\
{\
        if ( yychar >= 0 || ( yyr2[ yytmp ] >> 1 ) != 1 )\
        {\
                yyusercatd=catopen("yacc_user.cat", NL_CAT_LOCALE);\
                yyerror(catgets(yyusercatd,1,1,"syntax error - cannot backup" ));\
                YYERROR;\
        }\
        yychar = newtoken;\
        yystate = *yyps;\
        yylval = newvalue;\
        goto yynewstate;\
}
#endif
#define YYRECOVERING()  (!!yyerrflag)
#ifndef YYDEBUG
#       define YYDEBUG  1       /* make debugging available */
#endif

/*
** user known globals
*/
int yydebug;                    /* set to 1 to get debugging */

/*
** driver internal defines
*/
#define YYFLAG          (-1000)

#ifdef YYSPLIT
#   define YYSCODE { \
                        extern int (*_yyf[])(); \
                        register int yyret; \
                        if (_yyf[yytmp]) \
                            if ((yyret=(*_yyf[yytmp])()) == -2) \
                                    goto yyerrlab; \
                                else if (yyret>=0) return(yyret); \
                   }
#endif

/*
** global variables used by the parser
*/
YYSTYPE yyv[ YYMAXDEPTH ];      /* value stack */
int yys[ YYMAXDEPTH ];          /* state stack */

YYSTYPE *yypv;                  /* top of value stack */
YYSTYPE *yypvt;                 /* top of value stack for $vars */
int *yyps;                      /* top of state stack */

int yystate;                    /* current state */
int yytmp;                      /* extra var (lasts between blocks) */

int yynerrs;                    /* number of errors */
int yyerrflag;                  /* error recovery flag */
int yychar;                     /* current input token number */

#ifdef __cplusplus
 #ifdef _CPP_IOSTREAMS
  #include <iostream.h>
  extern void yyerror (char *); /* error message routine -- iostream version */
 #else
  #include <stdio.h>
  extern "C" void yyerror (char *); /* error message routine -- stdio version */
 #endif /* _CPP_IOSTREAMS */
 extern "C" int yylex(void);        /* return the next token */
#endif /* __cplusplus */


/*
** yyparse - return 0 if worked, 1 if syntax error not recovered from
*/
#ifdef __cplusplus
extern "C"
#endif /* __cplusplus */
int
yyparse()
{
        /*
        ** Initialize externals - yyparse may be called more than once
        */
        yypv = &yyv[-1];
        yyps = &yys[-1];
        yystate = 0;
        yytmp = 0;
        yynerrs = 0;
        yyerrflag = 0;
        yychar = -1;
#ifdef YACC_MSG
        yyusercatd=catopen("yacc_user.cat", NL_CAT_LOCALE);
#endif
        goto yystack;
        {
                register YYSTYPE *yy_pv;        /* top of value stack */
                register int *yy_ps;            /* top of state stack */
                register int yy_state;          /* current state */
                register int  yy_n;             /* internal state number info */

                /*
                ** get globals into registers.
                ** branch to here only if YYBACKUP was called.
                */
        yynewstate:
                yy_pv = yypv;
                yy_ps = yyps;
                yy_state = yystate;
                goto yy_newstate;

                /*
                ** get globals into registers.
                ** either we just started, or we just finished a reduction
                */
        yystack:
                yy_pv = yypv;
                yy_ps = yyps;
                yy_state = yystate;

                /*
                ** top of for (;;) loop while no reductions done
                */
        yy_stack:
                /*
                ** put a state and value onto the stacks
                */
#if YYDEBUG
                /*
                ** if debugging, look up token value in list of value vs.
                ** name pairs.  0 and negative (-1) are special values.
                ** Note: linear search is used since time is not a real
                ** consideration while debugging.
                */
                if ( yydebug )
                {
                        register int yy_i;

#if defined(__cplusplus) && defined(_CPP_IOSTREAMS)
                        cout << "State " << yy_state << " token ";
                        if ( yychar == 0 )
                                cout << "end-of-file" << endl;
                        else if ( yychar < 0 )
                                cout << "-none-" << endl;
#else
                        printf( "State %d, token ", yy_state );
                        if ( yychar == 0 )
                                printf( "end-of-file\n" );
                        else if ( yychar < 0 )
                                printf( "-none-\n" );
#endif /* defined(__cplusplus) && defined(_CPP_IOSTREAMS) */
                        else
                        {
                                for ( yy_i = 0; yytoks[yy_i].t_val >= 0;
                                        yy_i++ )
                                {
                                        if ( yytoks[yy_i].t_val == yychar )
                                                break;
                                }
#if defined(__cplusplus) && defined(_CPP_IOSTREAMS)
                                cout << yytoks[yy_i].t_name << endl;
#else
                                printf( "%s\n", yytoks[yy_i].t_name );
#endif /* defined(__cplusplus) && defined(_CPP_IOSTREAMS) */
                        }
                }
#endif /* YYDEBUG */
                if ( ++yy_ps >= &yys[ YYMAXDEPTH ] )    /* room on stack? */
                {
#ifndef YACC_MSG
                        yyerror( "yacc stack overflow" );
#else
                        yyerror(catgets(yyusercatd,1,2,"yacc stack overflow" ));
#endif
                        YYABORT;
                }
                *yy_ps = yy_state;
                *++yy_pv = yyval;

                /*
                ** we have a new state - find out what to do
                */
        yy_newstate:
                if ( ( yy_n = yypact[ yy_state ] ) <= YYFLAG )
                        goto yydefault;         /* simple state */
#if YYDEBUG
                /*
                ** if debugging, need to mark whether new token grabbed
                */
                yytmp = yychar < 0;
#endif
                if ( ( yychar < 0 ) && ( ( yychar = yylex() ) < 0 ) )
                        yychar = 0;             /* reached EOF */
#if YYDEBUG
                if ( yydebug && yytmp )
                {
                        register int yy_i;

#if defined(__cplusplus) && defined(_CPP_IOSTREAMS)
                        cout << "Received token " << endl;
                        if ( yychar == 0 )
                                cout << "end-of-file" << endl;
                        else if ( yychar < 0 )
                                cout << "-none-" << endl;
#else
                        printf( "Received token " );
                        if ( yychar == 0 )
                                printf( "end-of-file\n" );
                        else if ( yychar < 0 )
                                printf( "-none-\n" );
#endif /* defined(__cplusplus) && defined(_CPP_IOSTREAMS) */
                        else
                        {
                                for ( yy_i = 0; yytoks[yy_i].t_val >= 0;
                                        yy_i++ )
                                {
                                        if ( yytoks[yy_i].t_val == yychar )
                                                break;
                                }
#if defined(__cplusplus) && defined(_CPP_IOSTREAMS)
                                cout << yytoks[yy_i].t_name << endl;
#else
                                printf( "%s\n", yytoks[yy_i].t_name );
#endif /* defined(__cplusplus) && defined(_CPP_IOSTREAMS) */
                        }
                }
#endif /* YYDEBUG */
                if ( ( ( yy_n += yychar ) < 0 ) || ( yy_n >= YYLAST ) )
                        goto yydefault;
                if ( yychk[ yy_n = yyact[ yy_n ] ] == yychar )  /*valid shift*/
                {
                        yychar = -1;
                        yyval = yylval;
                        yy_state = yy_n;
                        if ( yyerrflag > 0 )
                                yyerrflag--;
                        goto yy_stack;
                }

        yydefault:
                if ( ( yy_n = yydef[ yy_state ] ) == -2 )
                {
#if YYDEBUG
                        yytmp = yychar < 0;
#endif
                        if ( ( yychar < 0 ) && ( ( yychar = yylex() ) < 0 ) )
                                yychar = 0;             /* reached EOF */
#if YYDEBUG
                        if ( yydebug && yytmp )
                        {
                                register int yy_i;

#if defined(__cplusplus) && defined(_CPP_IOSTREAMS)
                                cout << "Received token " << endl;
                                if ( yychar == 0 )
                                        cout << "end-of-file" << endl;
                                else if ( yychar < 0 )
                                        cout << "-none-" << endl;
#else
                                printf( "Received token " );
                                if ( yychar == 0 )
                                        printf( "end-of-file\n" );
                                else if ( yychar < 0 )
                                        printf( "-none-\n" );
#endif /* defined(__cplusplus) && defined(_CPP_IOSTREAMS) */
                                else
                                {
                                        for ( yy_i = 0;
                                                yytoks[yy_i].t_val >= 0;
                                                yy_i++ )
                                        {
                                                if ( yytoks[yy_i].t_val
                                                        == yychar )
                                                {
                                                        break;
                                                }
                                        }
#if defined(__cplusplus) && defined(_CPP_IOSTREAMS)
                                        cout << yytoks[yy_i].t_name << endl;
#else
                                        printf( "%s\n", yytoks[yy_i].t_name );
#endif /* defined(__cplusplus) && defined(_CPP_IOSTREAMS) */
                                }
                        }
#endif /* YYDEBUG */
                        /*
                        ** look through exception table
                        */
                        {
                                register int *yyxi = yyexca;

                                while ( ( *yyxi != -1 ) ||
                                        ( yyxi[1] != yy_state ) )
                                {
                                        yyxi += 2;
                                }
                                while ( ( *(yyxi += 2) >= 0 ) &&
                                        ( *yyxi != yychar ) )
                                        ;
                                if ( ( yy_n = yyxi[1] ) < 0 )
                                        YYACCEPT;
                        }
                }

                /*
                ** check for syntax error
                */
                if ( yy_n == 0 )        /* have an error */
                {
                        /* no worry about speed here! */
                        switch ( yyerrflag )
                        {
                        case 0:         /* new error */
#ifndef YACC_MSG
                                yyerror( "syntax error" );
#else
                                yyerror(catgets(yyusercatd,1,3,"syntax error" ));
#endif
                                goto skip_init;
                        yyerrlab:
                                /*
                                ** get globals into registers.
                                ** we have a user generated syntax type error
                                */
                                yy_pv = yypv;
                                yy_ps = yyps;
                                yy_state = yystate;
                                yynerrs++;
                        skip_init:
                        case 1:
                        case 2:         /* incompletely recovered error */
                                        /* try again... */
                                yyerrflag = 3;
                                /*
                                ** find state where "error" is a legal
                                ** shift action
                                */
                                while ( yy_ps >= yys )
                                {
                                        yy_n = yypact[ *yy_ps ] + YYERRCODE;
                                        if ( yy_n >= 0 && yy_n < YYLAST &&
                                                yychk[yyact[yy_n]] == YYERRCODE)                                        {
                                                /*
                                                ** simulate shift of "error"
                                                */
                                                yy_state = yyact[ yy_n ];
                                                goto yy_stack;
                                        }
                                        /*
                                        ** current state has no shift on
                                        ** "error", pop stack
                                        */
#if YYDEBUG
                                        if ( yydebug )
#if defined(__cplusplus) && defined(_CPP_IOSTREAMS)
                                            cout << "Error recovery pops state "
                                                 << (*yy_ps)
                                                 << ", uncovers state "
                                                 << yy_ps[-1] << endl;
#else
#       define _POP_ "Error recovery pops state %d, uncovers state %d\n"
                                                printf( _POP_, *yy_ps,
                                                        yy_ps[-1] );
#       undef _POP_
#endif /* defined(__cplusplus) && defined(_CPP_IOSTREAMS) */
#endif
                                        yy_ps--;
                                        yy_pv--;
                                }
                                /*
                                ** there is no state on stack with "error" as
                                ** a valid shift.  give up.
                                */
                                YYABORT;
                        case 3:         /* no shift yet; eat a token */
#if YYDEBUG
                                /*
                                ** if debugging, look up token in list of
                                ** pairs.  0 and negative shouldn't occur,
                                ** but since timing doesn't matter when
                                ** debugging, it doesn't hurt to leave the
                                ** tests here.
                                */
                                if ( yydebug )
                                {
                                        register int yy_i;

#if defined(__cplusplus) && defined(_CPP_IOSTREAMS)
                                        cout << "Error recovery discards ";
                                        if ( yychar == 0 )
                                            cout << "token end-of-file" << endl;
                                        else if ( yychar < 0 )
                                            cout << "token -none-" << endl;
#else
                                        printf( "Error recovery discards " );
                                        if ( yychar == 0 )
                                                printf( "token end-of-file\n" );
                                        else if ( yychar < 0 )
                                                printf( "token -none-\n" );
#endif /* defined(__cplusplus) && defined(_CPP_IOSTREAMS) */
                                        else
                                        {
                                                for ( yy_i = 0;
                                                        yytoks[yy_i].t_val >= 0;
                                                        yy_i++ )
                                                {
                                                        if ( yytoks[yy_i].t_val
                                                                == yychar )
                                                        {
                                                                break;
                                                        }
                                                }
#if defined(__cplusplus) && defined(_CPP_IOSTREAMS)
                                                cout << "token " <<
                                                    yytoks[yy_i].t_name <<
                                                    endl;
#else
                                                printf( "token %s\n",
                                                        yytoks[yy_i].t_name );
#endif /* defined(__cplusplus) && defined(_CPP_IOSTREAMS) */
                                        }
                                }
#endif /* YYDEBUG */
                                if ( yychar == 0 )      /* reached EOF. quit */
                                        YYABORT;
                                yychar = -1;
                                goto yy_newstate;
                        }
                }/* end if ( yy_n == 0 ) */
                /*
                ** reduction by production yy_n
                ** put stack tops, etc. so things right after switch
                */
#if YYDEBUG
                /*
                ** if debugging, print the string that is the user's
                ** specification of the reduction which is just about
                ** to be done.
                */
                if ( yydebug )
#if defined(__cplusplus) && defined(_CPP_IOSTREAMS)
                        cout << "Reduce by (" << yy_n << ") \"" <<
                            yyreds[ yy_n ] << "\"\n";
#else
                        printf( "Reduce by (%d) \"%s\"\n",
                                yy_n, yyreds[ yy_n ] );
#endif /* defined(__cplusplus) && defined(_CPP_IOSTREAMS) */
#endif
                yytmp = yy_n;                   /* value to switch over */
                yypvt = yy_pv;                  /* $vars top of value stack */
                /*
                ** Look in goto table for next state
                ** Sorry about using yy_state here as temporary
                ** register variable, but why not, if it works...
                ** If yyr2[ yy_n ] doesn't have the low order bit
                ** set, then there is no action to be done for
                ** this reduction.  So, no saving & unsaving of
                ** registers done.  The only difference between the
                ** code just after the if and the body of the if is
                ** the goto yy_stack in the body.  This way the test
                ** can be made before the choice of what to do is needed.
                */
                {
                        /* length of production doubled with extra bit */
                        register int yy_len = yyr2[ yy_n ];

                        if ( !( yy_len & 01 ) )
                        {
                                yy_len >>= 1;
                                yyval = ( yy_pv -= yy_len )[1]; /* $$ = $1 */
                                yy_state = yypgo[ yy_n = yyr1[ yy_n ] ] +
                                        *( yy_ps -= yy_len ) + 1;
                                if ( yy_state >= YYLAST ||
                                        yychk[ yy_state =
                                        yyact[ yy_state ] ] != -yy_n )
                                {
                                        yy_state = yyact[ yypgo[ yy_n ] ];
                                }
                                goto yy_stack;
                        }
                        yy_len >>= 1;
                        yyval = ( yy_pv -= yy_len )[1]; /* $$ = $1 */
                        yy_state = yypgo[ yy_n = yyr1[ yy_n ] ] +
                                *( yy_ps -= yy_len ) + 1;
                        if ( yy_state >= YYLAST ||
                                yychk[ yy_state = yyact[ yy_state ] ] != -yy_n )
                        {
                                yy_state = yyact[ yypgo[ yy_n ] ];
                        }
                }
                                        /* save until reenter driver code */
                yystate = yy_state;
                yyps = yy_ps;
                yypv = yy_pv;
        }
        /*
        ** code supplied by user is placed in this switch
        */

                switch(yytmp){

case 1:
# line 115 "kc.y"
{ 
	      SetupTableMan(); 
	    } /*NOTREACHED*/ break;
case 2:
# line 119 "kc.y"
{
	      if (verbose==1) {               /* Printing status report */
		printf("Species used (*=dyn. var.):\n");
		for(i=1; i<=NoOfSpec(); i++) {
		  GetSpecNo(i, name, &charge);
		  if (IsSpecInConstraint(name, charge)==0)
		    printf("  *%s(", name);
		  else
		    printf("   %s(", name);
		  if (charge==FLT_MAX)
		    printf(".");
		  else {
		    if (charge>0.0)
		      printf("%d+", (int)charge);
		    else if (charge<0.0)
		      printf("%d-", -(int)charge);
		  } /* else */
		  printf(")\n");
		} /* for i*/
		for(i=1; i<=NumOfDynVar(); i++) {
		  GetDynVarNo(i, name);
		  printf("  *%s\n", name);
		} /* for i */
		if (NumOfParameter()!=0) {
		  printf("Parameters declared:\n");
		  for(i=1; i<=NumOfParameter(); i++) {
		    GetParamNo(i, name, &charge, &side); /* abuse of side */
		    if (charge==0.0)
		      printf("  %s\n", name);
		    else {
		      if (side==1) 
			printf("  %s\n", name);
		      else {
			printf("  %s(", name);
			if (charge==FLT_MAX)
			  printf(".");
			else {
			  if (charge>0.0)
			    printf("%d+", (int)charge);
			  else if (charge<0.0)
			    printf("%d-", -(int)charge);
			} /* else */
			printf(")\n");
		      } /* else */
		    } /* else */
		  } /* for i */
		} /* if */
		printf("Constants declared:\n");
		for(i=1; i<=NumOfConstants(); i++) {
		  GetConstantNo(i, name);
		  printf("  %s = %e\n", name, GetConstant(name));
		} /* for i */
		for(i=1; i<=NumOfStrConst(); i++) {
		  GetStrConstNo(i, name, string);
		  printf("  %s = \"%s\"\n", name, string);
		} /* for i */
	      } /* if */
	      if (IsNonAutoSystem()==1) {
		NewDynVar("time");
		tmp=TreeCreate();
		TreeAssignConst(tmp, 1.0);
		NewExpr("time", tmp);
		TreeKill(tmp);
	      } /* if */
	    } /*NOTREACHED*/ break;
case 8:
# line 190 "kc.y"
{ if (yypvt[-0].strnum.flag==1)
		NewConstant(yypvt[-2].name, yypvt[-0].strnum.numb);
	      else {
		NewStrConst(yypvt[-2].name, yypvt[-0].strnum.name);
	      }
	    } /*NOTREACHED*/ break;
case 9:
# line 197 "kc.y"
{ temp=TreeEval(yypvt[-1].tree);
	      if (TreeGetError()==NoEval) 
		fprintf(stderr, "Unable to evaluate expression in line %d.\n", lineno);
	      else {
		yyval.strnum.flag=1;
		yyval.strnum.numb=temp;
	      }
	      TreeKill(yypvt[-1].tree);
	    } /*NOTREACHED*/ break;
case 10:
# line 207 "kc.y"
{ yyval.strnum.flag=2;
	      strcpy(yyval.strnum.name, yypvt[-2].name);
	    } /*NOTREACHED*/ break;
case 13:
# line 213 "kc.y"
{ NewPrintVar(yypvt[-0].name); } /*NOTREACHED*/ break;
case 14:
# line 215 "kc.y"
{ NewPrintConc(yypvt[-0].comp.name, yypvt[-0].comp.charge); } /*NOTREACHED*/ break;
case 16:
# line 218 "kc.y"
{ temp=TreeEval(yypvt[-10].tree);
	      if (TreeGetError()==NoEval)
		fprintf(stderr, "WARNING: Could not evaluate expr in line %d\n", lineno);
	      else {
		NewParameter(yypvt[-12].name, temp);
	      }
	      TreeKill(yypvt[-10].tree); 
	      temp=TreeEval(yypvt[-8].tree);
	      if (TreeGetError()==NoEval)
		fprintf(stderr, "WARNING: Could not evaluate expr in line %d\n", lineno);
	      else
		NewDeltaParam(yypvt[-12].name, temp);
	      TreeKill(yypvt[-8].tree);
	      temp1=TreeEval(yypvt[-6].tree);
	      if (TreeGetError()==NoEval)
		fprintf(stderr, "WARNING: Could not evaluate expr in line %d\n" , lineno);
	      TreeKill(yypvt[-6].tree);
	      temp2=TreeEval(yypvt[-2].tree);
	      if (TreeGetError()==NoEval)
		fprintf(stderr, "WARNING: Could not evaluate expr in line %d\n" , lineno);
	      TreeKill(yypvt[-2].tree);
	      NewLowHighPrefParam(yypvt[-12].name, temp, temp1, temp2);
	      temp=TreeEval(yypvt[-0].tree);
	      if (TreeGetError()==NoEval)
		fprintf(stderr, "WARNING: Could not evaluate expr in line %d\n", lineno);
	      TreeKill(yypvt[-0].tree);
	      NewDirectForParam(yypvt[-12].name, (int)temp);
	      temp=TreeEval(yypvt[-4].tree);
	      if (TreeGetError()==NoEval)
		fprintf(stderr, "WARNING: Could not evaluate expr in line %d\n", lineno);
	      TreeKill(yypvt[-4].tree);
	      NewDeltaParam(yypvt[-12].name, temp);
	    } /*NOTREACHED*/ break;
case 17:
# line 252 "kc.y"
{ temp=TreeEval(yypvt[-10].tree);
	      if (TreeGetError()==NoEval)
		fprintf(stderr, "WARNING: Could not evaluate expr in line %\n", lineno);
	      else
		NewParamConc(yypvt[-12].comp.name, yypvt[-12].comp.charge, temp);
	      temp=TreeEval(yypvt[-8].tree);
	      if (TreeGetError()==NoEval)
		fprintf(stderr, "WARNING: Could not evaluate expr in line %\n", lineno);
	      TreeKill(yypvt[-10].tree);
	      TreeKill(yypvt[-8].tree);
	      temp1=TreeEval(yypvt[-6].tree);
	      if (TreeGetError()==NoEval)
		fprintf(stderr, "WARNING: Could not evaluate expr in line %d\n" , lineno);
	      TreeKill(yypvt[-6].tree);
	      temp2=TreeEval(yypvt[-2].tree);
	      if (TreeGetError()==NoEval)
		fprintf(stderr, "WARNING: Could not evaluate expr in line %d\n" , lineno);
	      TreeKill(yypvt[-2].tree);
	      NewLowHighPrefConc(yypvt[-12].comp.name, yypvt[-12].comp.charge, temp, temp1, temp2);
	      temp=TreeEval(yypvt[-0].tree);
	      if (TreeGetError()==NoEval)
		fprintf(stderr, "WARNING: Could not evaluate expr in line %d\n", lineno);
	      TreeKill(yypvt[-4].tree);
	      NewDirectForConc(yypvt[-12].comp.name, yypvt[-12].comp.charge, (int)temp);
	      temp=TreeEval(yypvt[-4].tree);
	      if (TreeGetError()==NoEval)
		fprintf(stderr, "WARNING: Could not evaluate expr in line %d\n", lineno);
	      TreeKill(yypvt[-4].tree);
	      NewDeltaConc(yypvt[-12].comp.name, yypvt[-12].comp.charge, temp);
	    } /*NOTREACHED*/ break;
case 20:
# line 285 "kc.y"
{ NewReaction((int)value); 
	      side=1; 
	    } /*NOTREACHED*/ break;
case 21:
# line 289 "kc.y"
{ 
	      switch (yypvt[-0].flag) {
	      case '>': 
		AddReactionKind(GetCurrentReaction(), uni); 
		break;
	      case '<': 
		AddReactionKind(GetCurrentReaction(), bi); 
		break;
	      case '=': 
		AddReactionKind(GetCurrentReaction(), equi); 
		break;
	      }; /* switch */
	      side=-1;
	    } /*NOTREACHED*/ break;
case 23:
# line 305 "kc.y"
{ NewExpr(yypvt[-4].name, yypvt[-1].tree);
	      NewDynVar(yypvt[-4].name);
	      TreeKill(yypvt[-1].tree);
	    } /*NOTREACHED*/ break;
case 28:
# line 314 "kc.y"
{ if (strcmp(yypvt[-5].name, "c")!=0) 
		fprintf(stderr, "c expected!\n");
	      else {
		temp=TreeEval(yypvt[-0].tree);
		if (TreeGetError()==NoEval) 
		  fprintf(stderr, "Unable to evaluate expression in line %d.\n", lineno);
		else
		  NewPowerConst(GetCurrentReaction(), yypvt[-3].comp.name, yypvt[-3].comp.charge, temp, side); 
	      }; /* else */
	      TreeKill(yypvt[-0].tree);
	    } /*NOTREACHED*/ break;
case 29:
# line 326 "kc.y"
{ if (strcmp(yypvt[-2].name, "k")==0) {  
		NewRateConst(GetCurrentReaction(), -1, yypvt[-0].tree); 
	      } /* if */
	      else 
		if (strcmp(yypvt[-2].name, "v")==0) 
		  NewRateExpr(GetCurrentReaction(), -1, yypvt[-0].tree);
		else
		  fprintf(stderr, "Syntax error: v or k expected in line %d.\n", lineno);
	      TreeKill(yypvt[-0].tree);               
	    } /*NOTREACHED*/ break;
case 31:
# line 338 "kc.y"
{ if (strcmp(yypvt[-2].name, "K")!=0) 
		fprintf(stderr, "Syntax error: K expected in line %d.\n", lineno);
	      else {
		NewRateConst(GetCurrentReaction(), 0, yypvt[-0].tree);
	      }; /* else */
	      TreeKill(yypvt[-0].tree); 
	    } /*NOTREACHED*/ break;
case 32:
# line 345 "kc.y"
{ if (GetReactKind(GetCurrentReaction())==bi) (void) fprintf(stderr, "The reaction in line %d is two-ways.\n", lineno);
	    } /*NOTREACHED*/ break;
case 33:
# line 348 "kc.y"
{ if (GetReactKind(GetCurrentReaction())!=bi)   
		  (void) fprintf(stderr, "The reaction in line %d is a >one-way< reaction or >equilibrium.<\n", lineno);
		else {
		  if (strcmp(yypvt[-3].name, "k")==0)  
		    NewRateConst(GetCurrentReaction(), 1, yypvt[-1].tree);  
		  else
		    if (strcmp(yypvt[-3].name, "v")==0) 
		      NewRateExpr(GetCurrentReaction(), 1, yypvt[-1].tree);
		    else
		      fprintf(stderr, "Syntax error: v or k expected in line %d.\n", lineno);
	      }; /* else */
	      TreeKill(yypvt[-1].tree);
	    } /*NOTREACHED*/ break;
case 34:
# line 361 "kc.y"
{ yyval.flag='>'; } /*NOTREACHED*/ break;
case 35:
# line 362 "kc.y"
{ yyval.flag='<'; } /*NOTREACHED*/ break;
case 36:
# line 363 "kc.y"
{ yyval.flag='='; } /*NOTREACHED*/ break;
case 37:
# line 365 "kc.y"
{ SpecieInReaction(GetCurrentReaction(), yypvt[-0].comp.name, yypvt[-0].comp.charge);
	      NewCoeff(GetCurrentReaction(), yypvt[-0].comp.name, yypvt[-0].comp.charge, yypvt[-1].dval, side);
	      NewSpecie(yypvt[-0].comp.name, yypvt[-0].comp.charge);
	    } /*NOTREACHED*/ break;
case 38:
# line 370 "kc.y"
{ SpecieInReaction(GetCurrentReaction(), yypvt[-0].comp.name, yypvt[-0].comp.charge);
	      NewCoeff(GetCurrentReaction(), yypvt[-0].comp.name, yypvt[-0].comp.charge, yypvt[-1].dval, side);
	      NewSpecie(yypvt[-0].comp.name, yypvt[-0].comp.charge);
	    } /*NOTREACHED*/ break;
case 39:
# line 375 "kc.y"
{ 
	      (void) strcpy(yyval.comp.name, yypvt[-1].name);
	      yyval.comp.charge=yypvt[-0].dval;
	      yyval.comp.concs=no;
	    } /*NOTREACHED*/ break;
case 40:
# line 381 "kc.y"
{ yyval.dval=yypvt[-1].dval; } /*NOTREACHED*/ break;
case 41:
# line 382 "kc.y"
{ 
		yyval.dval=0.0;
	      } /*NOTREACHED*/ break;
case 42:
# line 386 "kc.y"
{ charge=FLT_MAX;
	      yyval.dval=FLT_MAX;
	    } /*NOTREACHED*/ break;
case 43:
# line 390 "kc.y"
{ if (yypvt[-1].oper=='+') {
		yyval.dval=yypvt[-0].dval;
		charge=yypvt[-0].dval;
	      } /* if */
	      else {
		charge=-yypvt[-0].dval;
		yyval.dval=-yypvt[-0].dval;
	      }; /* else */
	    } /*NOTREACHED*/ break;
case 44:
# line 400 "kc.y"
{ if (yypvt[-0].oper=='+')
		yyval.dval=1.0;
	      else
		yyval.dval=-1.0;
	    } /*NOTREACHED*/ break;
case 45:
# line 405 "kc.y"
{ yyval.oper = '+'; } /*NOTREACHED*/ break;
case 46:
# line 406 "kc.y"
{ yyval.oper = '-'; } /*NOTREACHED*/ break;
case 47:
# line 408 "kc.y"
{ coeff=value;
	      yyval.dval=yypvt[-0].dval;  
	    } /*NOTREACHED*/ break;
case 48:
# line 412 "kc.y"
{ coeff=1.0;
	      yyval.dval=1.0;
	    } /*NOTREACHED*/ break;
case 51:
# line 418 "kc.y"
{ temp=TreeEval(yypvt[-1].tree); 
	      switch (yypvt[-3].flag) {
	      case 1: 
		if (TreeGetError()==NoEval) 
		  fprintf(stderr, "Unable to evaluate expression.\n");
		else  
		  NewBeginConc(yypvt[-4].comp.name, yypvt[-4].comp.charge, temp); 
		break;
	      case 0: 
		NewConstraint(yypvt[-4].comp.name, yypvt[-4].comp.charge, yypvt[-1].tree);
		break;
	      };
	      TreeKill(yypvt[-1].tree);
	    } /*NOTREACHED*/ break;
case 52:
# line 433 "kc.y"
{ temp=TreeEval(yypvt[-1].tree);
	      NewSpecConst(yypvt[-4].comp.name, yypvt[-4].comp.charge, yypvt[-6].name, temp);
	      if (GetError()==NotFound) 
		if (yypvt[-4].comp.charge==0.0) 
		  NewDynVarConst(yypvt[-4].comp.name, yypvt[-6].name, temp);
	      TreeKill(yypvt[-1].tree);
	    } /*NOTREACHED*/ break;
case 53:
# line 441 "kc.y"
{ temp=TreeEval(yypvt[-1].tree);
	      if (TreeGetError()==NoEval)
		yyerror("Unable to evaluate expression");
	      else {
		if (yypvt[-3].flag==1) {
		  NewInitValue(yypvt[-4].name, temp);
		  TreeKill(yypvt[-1].tree);
		} else
		  yyerror("(0) expected");
	      }
	    } /*NOTREACHED*/ break;
case 54:
# line 453 "kc.y"
{ yyval.tree=TreeCreate();
	      TreeCpy(yyval.tree, yypvt[-0].tree);
	      TreeSign(yyval.tree);
	      TreeKill(yypvt[-0].tree);
	    } /*NOTREACHED*/ break;
case 55:
# line 459 "kc.y"
{ yyval.tree=TreeCreate(); 
	      TreeAdd(yypvt[-2].tree, yypvt[-0].tree);
	      TreeCpy(yyval.tree, yypvt[-2].tree);
	      TreeKill(yypvt[-0].tree);
	      TreeKill(yypvt[-2].tree);
	    } /*NOTREACHED*/ break;
case 56:
# line 466 "kc.y"
{ yyval.tree=TreeCreate(); 
	      TreeSub(yypvt[-2].tree, yypvt[-0].tree);
	      TreeCpy(yyval.tree, yypvt[-2].tree);
	      TreeKill(yypvt[-0].tree);
	      TreeKill(yypvt[-2].tree);
	    } /*NOTREACHED*/ break;
case 57:
# line 473 "kc.y"
{
	      yyval.tree=TreeCreate();
	      TreeMul(yypvt[-2].tree, yypvt[-0].tree);
	      TreeCpy(yyval.tree, yypvt[-2].tree);
	      TreeKill(yypvt[-2].tree);
	      TreeKill(yypvt[-0].tree);
	    } /*NOTREACHED*/ break;
case 58:
# line 481 "kc.y"
{
	      yyval.tree=TreeCreate();
	      TreeDiv(yypvt[-2].tree, yypvt[-0].tree);
	      TreeCpy(yyval.tree, yypvt[-2].tree);
	      TreeKill(yypvt[-2].tree);
	      TreeKill(yypvt[-0].tree);
	    } /*NOTREACHED*/ break;
case 59:
# line 489 "kc.y"
{ yyval.tree=TreeCreate();
	      TreeCpy(yyval.tree, yypvt[-2].tree);
	      TreePow(yyval.tree, yypvt[-0].tree);
	      TreeKill(yypvt[-2].tree);
	      TreeKill(yypvt[-0].tree);
	    } /*NOTREACHED*/ break;
case 60:
# line 496 "kc.y"
{ yyval.tree=TreeCreate();
	      TreeCpy(yyval.tree, yypvt[-1].tree); 
	      TreeKill(yypvt[-1].tree);
	    } /*NOTREACHED*/ break;
case 61:
# line 501 "kc.y"
{
	     yyval.tree=TreeCreate();
	     temp=GetConstant(yypvt[-0].name); 
	     if (GetError()==NotFound) {
	       TreeAssignVar(yyval.tree, yypvt[-0].name); 
	       if (strcmp(yypvt[-0].name, "time")==0)
		 NonAutoSystem();
	     } else
	       TreeAssignConst(yyval.tree, temp);
	    } /*NOTREACHED*/ break;
case 62:
# line 512 "kc.y"
{ yyval.tree=TreeCreate();
	      TreeAssignConst(yyval.tree, yypvt[-0].dval);
	    } /*NOTREACHED*/ break;
case 63:
# line 516 "kc.y"
{  
	      yyval.tree=TreeCreate(); 
	      if (yypvt[-0].flag==1) {
		temp=GetBeginConc(yypvt[-1].comp.name, yypvt[-1].comp.charge);
		if (GetError()==NoError)
		  TreeAssignConst(yyval.tree, temp);
		else
		  fprintf(stderr, "[%s(%e)] not found.\n", yypvt[-1].comp.name, yypvt[-1].comp.charge);
		flag='1';
	      } /* if */
	      else { 
		flag='0';
		RenameSpec(name, yypvt[-1].comp.name, yypvt[-1].comp.charge);
		TreeAssignVar(yyval.tree, name);
	      }; /* else */
	    } /*NOTREACHED*/ break;
case 64:
# line 533 "kc.y"
{ yyval.tree=TreeCreate();
	      TreeCpy(yyval.tree, yypvt[-1].tree);
	      TreeApplyFunc(&yyval.tree, yypvt[-3].func);
	      TreeKill(yypvt[-1].tree);
	    } /*NOTREACHED*/ break;
case 65:
# line 538 "kc.y"
{ yyval.func=Exp; } /*NOTREACHED*/ break;
case 66:
# line 539 "kc.y"
{ yyval.func=Ln; } /*NOTREACHED*/ break;
case 67:
# line 540 "kc.y"
{ yyval.func=Log; } /*NOTREACHED*/ break;
case 68:
# line 541 "kc.y"
{ yyval.func=Sin; } /*NOTREACHED*/ break;
case 69:
# line 542 "kc.y"
{ yyval.func=Cos; } /*NOTREACHED*/ break;
case 70:
# line 543 "kc.y"
{ yyval.func=Tan; } /*NOTREACHED*/ break;
case 71:
# line 544 "kc.y"
{ yyval.func=Sinh; } /*NOTREACHED*/ break;
case 72:
# line 545 "kc.y"
{ yyval.func=Cosh; } /*NOTREACHED*/ break;
case 73:
# line 546 "kc.y"
{ yyval.func=Tanh; } /*NOTREACHED*/ break;
case 74:
# line 547 "kc.y"
{ yyval.func=Asin; } /*NOTREACHED*/ break;
case 75:
# line 548 "kc.y"
{ yyval.func=Acos; } /*NOTREACHED*/ break;
case 76:
# line 549 "kc.y"
{ yyval.func=Atan; } /*NOTREACHED*/ break;
case 77:
# line 550 "kc.y"
{ yyval.func=Asinh; } /*NOTREACHED*/ break;
case 78:
# line 551 "kc.y"
{ yyval.func=Acosh; } /*NOTREACHED*/ break;
case 79:
# line 552 "kc.y"
{ yyval.func=Atanh; } /*NOTREACHED*/ break;
case 80:
# line 554 "kc.y"
{ 
	      (void) strcpy(yyval.comp.name, yypvt[-1].comp.name);
	      yyval.comp.charge=yypvt[-1].comp.charge;
	      yyval.comp.concs=ord;
	    } /*NOTREACHED*/ break;
case 81:
# line 559 "kc.y"
{ yyval.flag=0; } /*NOTREACHED*/ break;
case 82:
# line 561 "kc.y"
{ yyval.flag=1; } /*NOTREACHED*/ break;
}


        goto yystack;           /* reset registers in driver code */
}

# line 562 "kc.y"

#include "lex.c"
