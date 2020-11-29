/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton implementation for Bison's Yacc-like parsers in C

   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005, 2006
   Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA 02110-1301, USA.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "2.3"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Using locations.  */
#define YYLSP_NEEDED 0



/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     names = 258,
     leftarr = 259,
     rightarr = 260,
     numbers = 261,
     param = 262,
     print = 263,
     powop = 264,
     semicolon = 265,
     quotation = 266,
     equal = 267,
     colon = 268,
     R = 269,
     E = 270,
     powc = 271,
     leftpar = 272,
     rightpar = 273,
     oneway = 274,
     twoways = 275,
     plus = 276,
     minus = 277,
     multi = 278,
     pdiv = 279,
     comma = 280,
     leftconc = 281,
     rightconc = 282,
     time0 = 283,
     K = 284,
     radical = 285,
     V = 286,
     prime = 287,
     fun_exp = 288,
     fun_log = 289,
     fun_ln = 290,
     fun_sin = 291,
     fun_cos = 292,
     fun_tan = 293,
     fun_sinh = 294,
     fun_cosh = 295,
     fun_tanh = 296,
     fun_asin = 297,
     fun_acos = 298,
     fun_atan = 299,
     fun_asinh = 300,
     fun_acosh = 301,
     fun_atanh = 302,
     UMINUS = 303
   };
#endif
/* Tokens.  */
#define names 258
#define leftarr 259
#define rightarr 260
#define numbers 261
#define param 262
#define print 263
#define powop 264
#define semicolon 265
#define quotation 266
#define equal 267
#define colon 268
#define R 269
#define E 270
#define powc 271
#define leftpar 272
#define rightpar 273
#define oneway 274
#define twoways 275
#define plus 276
#define minus 277
#define multi 278
#define pdiv 279
#define comma 280
#define leftconc 281
#define rightconc 282
#define time0 283
#define K 284
#define radical 285
#define V 286
#define prime 287
#define fun_exp 288
#define fun_log 289
#define fun_ln 290
#define fun_sin 291
#define fun_cos 292
#define fun_tan 293
#define fun_sinh 294
#define fun_cosh 295
#define fun_tanh 296
#define fun_asin 297
#define fun_acos 298
#define fun_atan 299
#define fun_asinh 300
#define fun_acosh 301
#define fun_atanh 302
#define UMINUS 303




/* Copy the first part of user declarations.  */
#line 1 "kc.y"
 
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


/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* Enabling the token table.  */
#ifndef YYTOKEN_TABLE
# define YYTOKEN_TABLE 0
#endif

#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
#line 47 "kc.y"
{
  double      dval;
  char        oper;
  char        name[STRING_LENGTH];
  compound    comp;
  char        flag;
  Tree        tree;
  Function    func;
  Strnum      strnum;
}
/* Line 193 of yacc.c.  */
#line 249 "y.tab.c"
	YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif



/* Copy the second part of user declarations.  */


/* Line 216 of yacc.c.  */
#line 262 "y.tab.c"

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#elif (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
typedef signed char yytype_int8;
#else
typedef short int yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(msgid) dgettext ("bison-runtime", msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(msgid) msgid
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(e) ((void) (e))
#else
# define YYUSE(e) /* empty */
#endif

/* Identity function, used to suppress warnings about constant conditions.  */
#ifndef lint
# define YYID(n) (n)
#else
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static int
YYID (int i)
#else
static int
YYID (i)
    int i;
#endif
{
  return i;
}
#endif

#if ! defined yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#     ifndef _STDLIB_H
#      define _STDLIB_H 1
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (YYID (0))
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined _STDLIB_H \
       && ! ((defined YYMALLOC || defined malloc) \
	     && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef _STDLIB_H
#    define _STDLIB_H 1
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
	 || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss;
  YYSTYPE yyvs;
  };

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(To, From, Count) \
      __builtin_memcpy (To, From, (Count) * sizeof (*(From)))
#  else
#   define YYCOPY(To, From, Count)		\
      do					\
	{					\
	  YYSIZE_T yyi;				\
	  for (yyi = 0; yyi < (Count); yyi++)	\
	    (To)[yyi] = (From)[yyi];		\
	}					\
      while (YYID (0))
#  endif
# endif

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack)					\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack, Stack, yysize);				\
	Stack = &yyptr->Stack;						\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (YYID (0))

#endif

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  3
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   301

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  49
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  34
/* YYNRULES -- Number of rules.  */
#define YYNRULES  83
/* YYNRULES -- Number of states.  */
#define YYNSTATES  175

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   303

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      45,    46,    47,    48
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const yytype_uint16 yyprhs[] =
{
       0,     0,     3,     4,     9,    12,    13,    15,    17,    21,
      25,    28,    33,    37,    39,    41,    43,    47,    61,    75,
      77,    80,    81,    82,    92,    98,   102,   104,   108,   110,
     117,   118,   125,   129,   130,   135,   137,   139,   141,   144,
     149,   152,   156,   157,   159,   162,   164,   166,   168,   170,
     171,   173,   176,   182,   190,   196,   199,   203,   207,   211,
     215,   219,   223,   225,   227,   230,   235,   237,   239,   241,
     243,   245,   247,   249,   251,   253,   255,   257,   259,   261,
     263,   265,   269,   270
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int8 yyrhs[] =
{
      50,     0,    -1,    -1,    51,    52,    60,    77,    -1,    52,
      53,    -1,    -1,    54,    -1,    58,    -1,     8,    56,    10,
      -1,     3,    12,    55,    -1,    79,    10,    -1,    11,     3,
      11,    10,    -1,    56,    25,    57,    -1,    57,    -1,     3,
      -1,    81,    -1,     7,    59,    10,    -1,     3,    12,    79,
      25,    79,    25,    79,    25,    79,    25,    79,    25,    79,
      -1,    81,    12,    79,    25,    79,    25,    79,    25,    79,
      25,    79,    25,    79,    -1,    61,    -1,    60,    61,    -1,
      -1,    -1,     6,    62,    13,    71,    70,    63,    71,    10,
      64,    -1,     3,    32,    12,    79,    10,    -1,    65,    10,
      67,    -1,    67,    -1,    65,    10,    66,    -1,    66,    -1,
       3,    17,    72,    18,    12,    79,    -1,    -1,     4,    12,
      79,    68,    10,    69,    -1,     3,    12,    79,    -1,    -1,
       5,    12,    79,    10,    -1,    19,    -1,    20,    -1,    12,
      -1,    76,    72,    -1,    71,    21,    76,    72,    -1,     3,
      73,    -1,    17,    74,    18,    -1,    -1,    30,    -1,    75,
       6,    -1,    75,    -1,    21,    -1,    22,    -1,     6,    -1,
      -1,    78,    -1,    77,    78,    -1,    81,    82,    12,    79,
      10,    -1,     3,    17,    72,    18,    12,    79,    10,    -1,
       3,    82,    12,    79,    10,    -1,    22,    79,    -1,    79,
      21,    79,    -1,    79,    22,    79,    -1,    79,    23,    79,
      -1,    79,    24,    79,    -1,    79,     9,    79,    -1,    17,
      79,    18,    -1,     3,    -1,     6,    -1,    81,    82,    -1,
      80,    17,    79,    18,    -1,    33,    -1,    35,    -1,    34,
      -1,    36,    -1,    37,    -1,    38,    -1,    39,    -1,    40,
      -1,    41,    -1,    42,    -1,    43,    -1,    44,    -1,    45,
      -1,    46,    -1,    47,    -1,    26,    72,    27,    -1,    -1,
      28,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   116,   116,   116,   185,   186,   187,   188,   189,   190,
     197,   207,   211,   212,   213,   215,   217,   218,   252,   283,
     284,   286,   290,   285,   305,   310,   311,   312,   313,   314,
     327,   326,   338,   346,   348,   362,   363,   364,   365,   370,
     375,   381,   383,   386,   390,   400,   406,   407,   408,   413,
     416,   417,   418,   433,   441,   453,   459,   466,   473,   481,
     489,   496,   501,   512,   516,   533,   539,   540,   541,   542,
     543,   544,   545,   546,   547,   548,   549,   550,   551,   552,
     553,   554,   560,   561
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "names", "leftarr", "rightarr",
  "numbers", "param", "print", "powop", "semicolon", "quotation", "equal",
  "colon", "R", "E", "powc", "leftpar", "rightpar", "oneway", "twoways",
  "plus", "minus", "multi", "pdiv", "comma", "leftconc", "rightconc",
  "time0", "K", "radical", "V", "prime", "fun_exp", "fun_log", "fun_ln",
  "fun_sin", "fun_cos", "fun_tan", "fun_sinh", "fun_cosh", "fun_tanh",
  "fun_asin", "fun_acos", "fun_atan", "fun_asinh", "fun_acosh",
  "fun_atanh", "UMINUS", "$accept", "system", "@1", "declars", "declar",
  "constant", "strnumb", "printlist", "prn_entry", "parameter", "pentry",
  "reactions", "reaction", "@2", "@3", "reackonst", "powconsts",
  "powconst", "rateconsts", "@4", "ropt", "kind", "substs", "subst",
  "charge", "size", "sign", "coeff", "constants", "const", "expr",
  "function", "conc", "timeopt", 0
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,   290,   291,   292,   293,   294,
     295,   296,   297,   298,   299,   300,   301,   302,   303
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    49,    51,    50,    52,    52,    53,    53,    53,    54,
      55,    55,    56,    56,    57,    57,    58,    59,    59,    60,
      60,    62,    63,    61,    61,    64,    64,    65,    65,    66,
      68,    67,    67,    69,    69,    70,    70,    70,    71,    71,
      72,    73,    73,    74,    74,    74,    75,    75,    76,    76,
      77,    77,    78,    78,    78,    79,    79,    79,    79,    79,
      79,    79,    79,    79,    79,    79,    80,    80,    80,    80,
      80,    80,    80,    80,    80,    80,    80,    80,    80,    80,
      80,    81,    82,    82
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     0,     4,     2,     0,     1,     1,     3,     3,
       2,     4,     3,     1,     1,     1,     3,    13,    13,     1,
       2,     0,     0,     9,     5,     3,     1,     3,     1,     6,
       0,     6,     3,     0,     4,     1,     1,     1,     2,     4,
       2,     3,     0,     1,     2,     1,     1,     1,     1,     0,
       1,     2,     5,     7,     5,     2,     3,     3,     3,     3,
       3,     3,     1,     1,     2,     4,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     3,     0,     1
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       2,     0,     5,     1,     0,     0,    21,     0,     0,     4,
       6,     7,     0,    19,     0,     0,     0,     0,     0,     0,
       0,    14,     0,    13,    15,    82,    20,     3,    50,    82,
      62,    63,     0,     0,     0,    66,    68,    67,    69,    70,
      71,    72,    73,    74,    75,    76,    77,    78,    79,    80,
       9,     0,     0,    82,     0,    49,     0,    42,     0,    16,
       0,     8,     0,     0,    83,     0,    82,    51,     0,     0,
       0,    55,     0,    10,     0,     0,     0,     0,     0,    64,
       0,    48,     0,     0,     0,     0,    40,    81,     0,    12,
       0,     0,     0,     0,    61,    60,    56,    57,    58,    59,
       0,    24,    37,    35,    36,    49,    22,    38,     0,    46,
      47,    43,     0,    45,     0,     0,     0,     0,    11,    65,
       0,    49,     0,    41,    44,     0,     0,    54,    52,    39,
       0,     0,     0,     0,     0,     0,     0,    53,     0,     0,
      23,     0,    28,    26,     0,     0,     0,     0,     0,     0,
       0,     0,    32,     0,    30,    27,    25,     0,     0,     0,
       0,     0,     0,     0,    33,     0,     0,    29,     0,    31,
      17,    18,     0,     0,    34
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,     1,     2,     4,     9,    10,    50,    22,    23,    11,
      19,    12,    13,    16,   121,   140,   141,   142,   143,   160,
     169,   106,    82,    58,    86,   112,   113,    83,    27,    28,
      51,    52,    53,    65
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -60
static const yytype_int16 yypact[] =
{
     -60,    30,   -60,   -60,    12,    -6,   -60,     6,     7,   -60,
     -60,   -60,     5,   -60,   120,    24,    55,    26,    34,    59,
      60,   -60,     3,   -60,   -60,    45,   -60,     9,   -60,    46,
     -60,   -60,    77,   165,   165,   -60,   -60,   -60,   -60,   -60,
     -60,   -60,   -60,   -60,   -60,   -60,   -60,   -60,   -60,   -60,
     -60,   214,    75,    46,   165,    90,   165,    97,    70,   -60,
     165,   -60,     7,    34,   -60,   104,    -3,   -60,   105,   107,
     230,   -60,   165,   -60,   165,   165,   165,   165,   165,   -60,
     234,   -60,   108,    34,    25,    57,   -60,   -60,    31,   -60,
     101,   165,   165,   111,   -60,   -60,    67,    67,   113,   113,
     241,   -60,   -60,   -60,   -60,    90,   -60,   -60,   165,   -60,
     -60,   -60,   116,   130,   165,   128,   251,   257,   -60,   -60,
      34,    90,    42,   -60,   -60,    61,   165,   -60,   -60,   -60,
      50,   165,   165,   261,    13,    80,    85,   -60,    83,   129,
     -60,   133,   -60,   -60,   165,   165,   165,    34,   165,    13,
     126,   171,   277,   127,   277,   -60,   -60,   165,   165,   132,
     142,   204,   209,   165,   164,   165,   165,   277,   158,   -60,
     277,   277,   165,   267,   -60
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
     -60,   -60,   -60,   -60,   -60,   -60,   -60,   -60,    76,   -60,
     -60,   -60,   160,   -60,   -60,   -60,   -60,    27,    28,   -60,
     -60,   -60,    52,   -59,   -60,   -60,   -60,    69,   -60,   148,
     -33,   -60,    -5,   -24
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -1
static const yytype_uint8 yytable[] =
{
      70,    71,    20,    24,    90,    68,    14,    29,    25,    17,
      21,     6,    66,    61,    63,     5,   138,   139,     6,     7,
       8,    80,    29,    84,   107,    64,    15,    88,    62,    79,
       3,    18,    18,    18,    72,    18,    54,    57,    56,    95,
      72,    96,    97,    98,    99,   100,    74,    75,    76,    77,
     108,    72,    74,    75,    76,    77,   114,    24,   116,   117,
     134,   129,    63,    74,    75,    76,    77,   131,    55,    59,
      72,   105,    60,    64,    64,   122,    72,    15,   109,   110,
      69,   125,    74,    75,    76,    77,   132,   111,   153,    72,
      76,    77,    78,   133,    72,   146,    81,    87,   135,   136,
     147,    74,    75,    76,    77,   144,    74,    75,    76,    77,
     145,   150,   151,   152,    85,   154,    91,    92,    93,   115,
     102,   118,    72,    30,   161,   162,    31,   103,   104,   105,
     167,    32,   170,   171,   123,    72,   124,    33,    89,   173,
     126,   148,    34,   149,   163,   159,    18,    74,    75,    76,
      77,   157,   164,    35,    36,    37,    38,    39,    40,    41,
      42,    43,    44,    45,    46,    47,    48,    49,    30,   168,
     172,    31,    26,   130,   120,    67,   155,   156,     0,     0,
      72,     0,    33,     0,     0,     0,     0,    34,     0,     0,
       0,    18,    74,    75,    76,    77,   158,     0,    35,    36,
      37,    38,    39,    40,    41,    42,    43,    44,    45,    46,
      47,    48,    49,    72,     0,     0,     0,     0,    72,     0,
       0,     0,     0,    72,    73,    74,    75,    76,    77,   165,
      74,    75,    76,    77,   166,    74,    75,    76,    77,    72,
       0,     0,     0,    72,   101,     0,     0,     0,    94,     0,
      72,    74,    75,    76,    77,    74,    75,    76,    77,   119,
      72,   127,    74,    75,    76,    77,    72,   128,     0,     0,
      72,   137,    74,    75,    76,    77,    72,   174,    74,    75,
      76,    77,    74,    75,    76,    77,    72,     0,    74,    75,
      76,    77,     0,     0,     0,     0,     0,     0,    74,    75,
      76,    77
};

static const yytype_int16 yycheck[] =
{
      33,    34,     7,     8,    63,    29,    12,    12,     3,     3,
       3,     6,     3,    10,    17,     3,     3,     4,     6,     7,
       8,    54,    27,    56,    83,    28,    32,    60,    25,    53,
       0,    26,    26,    26,     9,    26,    12,     3,    12,    72,
       9,    74,    75,    76,    77,    78,    21,    22,    23,    24,
      25,     9,    21,    22,    23,    24,    25,    62,    91,    92,
      10,   120,    17,    21,    22,    23,    24,    25,    13,    10,
       9,    21,    12,    28,    28,   108,     9,    32,    21,    22,
       3,   114,    21,    22,    23,    24,    25,    30,   147,     9,
      23,    24,    17,   126,     9,    12,     6,    27,   131,   132,
      17,    21,    22,    23,    24,    25,    21,    22,    23,    24,
      25,   144,   145,   146,    17,   148,    12,    12,    11,    18,
      12,    10,     9,     3,   157,   158,     6,    19,    20,    21,
     163,    11,   165,   166,    18,     9,     6,    17,    62,   172,
      12,    12,    22,    10,    12,    18,    26,    21,    22,    23,
      24,    25,    10,    33,    34,    35,    36,    37,    38,    39,
      40,    41,    42,    43,    44,    45,    46,    47,     3,     5,
      12,     6,    12,   121,   105,    27,   149,   149,    -1,    -1,
       9,    -1,    17,    -1,    -1,    -1,    -1,    22,    -1,    -1,
      -1,    26,    21,    22,    23,    24,    25,    -1,    33,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      45,    46,    47,     9,    -1,    -1,    -1,    -1,     9,    -1,
      -1,    -1,    -1,     9,    10,    21,    22,    23,    24,    25,
      21,    22,    23,    24,    25,    21,    22,    23,    24,     9,
      -1,    -1,    -1,     9,    10,    -1,    -1,    -1,    18,    -1,
       9,    21,    22,    23,    24,    21,    22,    23,    24,    18,
       9,    10,    21,    22,    23,    24,     9,    10,    -1,    -1,
       9,    10,    21,    22,    23,    24,     9,    10,    21,    22,
      23,    24,    21,    22,    23,    24,     9,    -1,    21,    22,
      23,    24,    -1,    -1,    -1,    -1,    -1,    -1,    21,    22,
      23,    24
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,    50,    51,     0,    52,     3,     6,     7,     8,    53,
      54,    58,    60,    61,    12,    32,    62,     3,    26,    59,
      81,     3,    56,    57,    81,     3,    61,    77,    78,    81,
       3,     6,    11,    17,    22,    33,    34,    35,    36,    37,
      38,    39,    40,    41,    42,    43,    44,    45,    46,    47,
      55,    79,    80,    81,    12,    13,    12,     3,    72,    10,
      12,    10,    25,    17,    28,    82,     3,    78,    82,     3,
      79,    79,     9,    10,    21,    22,    23,    24,    17,    82,
      79,     6,    71,    76,    79,    17,    73,    27,    79,    57,
      72,    12,    12,    11,    18,    79,    79,    79,    79,    79,
      79,    10,    12,    19,    20,    21,    70,    72,    25,    21,
      22,    30,    74,    75,    25,    18,    79,    79,    10,    18,
      76,    63,    79,    18,     6,    79,    12,    10,    10,    72,
      71,    25,    25,    79,    10,    79,    79,    10,     3,     4,
      64,    65,    66,    67,    25,    25,    12,    17,    12,    10,
      79,    79,    79,    72,    79,    66,    67,    25,    25,    18,
      68,    79,    79,    12,    10,    25,    25,    79,     5,    69,
      79,    79,    12,    79,    10
};

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(yychar = YYEMPTY)
#define YYEMPTY		(-2)
#define YYEOF		0

#define YYACCEPT	goto yyacceptlab
#define YYABORT		goto yyabortlab
#define YYERROR		goto yyerrorlab


/* Like YYERROR except do call yyerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */

#define YYFAIL		goto yyerrlab

#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)					\
do								\
  if (yychar == YYEMPTY && yylen == 1)				\
    {								\
      yychar = (Token);						\
      yylval = (Value);						\
      yytoken = YYTRANSLATE (yychar);				\
      YYPOPSTACK (1);						\
      goto yybackup;						\
    }								\
  else								\
    {								\
      yyerror (YY_("syntax error: cannot back up")); \
      YYERROR;							\
    }								\
while (YYID (0))


#define YYTERROR	1
#define YYERRCODE	256


/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

#define YYRHSLOC(Rhs, K) ((Rhs)[K])
#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)				\
    do									\
      if (YYID (N))                                                    \
	{								\
	  (Current).first_line   = YYRHSLOC (Rhs, 1).first_line;	\
	  (Current).first_column = YYRHSLOC (Rhs, 1).first_column;	\
	  (Current).last_line    = YYRHSLOC (Rhs, N).last_line;		\
	  (Current).last_column  = YYRHSLOC (Rhs, N).last_column;	\
	}								\
      else								\
	{								\
	  (Current).first_line   = (Current).last_line   =		\
	    YYRHSLOC (Rhs, 0).last_line;				\
	  (Current).first_column = (Current).last_column =		\
	    YYRHSLOC (Rhs, 0).last_column;				\
	}								\
    while (YYID (0))
#endif


/* YY_LOCATION_PRINT -- Print the location on the stream.
   This macro was not mandated originally: define only if we know
   we won't break user code: when these are the locations we know.  */

#ifndef YY_LOCATION_PRINT
# if defined YYLTYPE_IS_TRIVIAL && YYLTYPE_IS_TRIVIAL
#  define YY_LOCATION_PRINT(File, Loc)			\
     fprintf (File, "%d.%d-%d.%d",			\
	      (Loc).first_line, (Loc).first_column,	\
	      (Loc).last_line,  (Loc).last_column)
# else
#  define YY_LOCATION_PRINT(File, Loc) ((void) 0)
# endif
#endif


/* YYLEX -- calling `yylex' with the right arguments.  */

#ifdef YYLEX_PARAM
# define YYLEX yylex (YYLEX_PARAM)
#else
# define YYLEX yylex ()
#endif

/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)			\
do {						\
  if (yydebug)					\
    YYFPRINTF Args;				\
} while (YYID (0))

# define YY_SYMBOL_PRINT(Title, Type, Value, Location)			  \
do {									  \
  if (yydebug)								  \
    {									  \
      YYFPRINTF (stderr, "%s ", Title);					  \
      yy_symbol_print (stderr,						  \
		  Type, Value); \
      YYFPRINTF (stderr, "\n");						  \
    }									  \
} while (YYID (0))


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_value_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# else
  YYUSE (yyoutput);
# endif
  switch (yytype)
    {
      default:
	break;
    }
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (yytype < YYNTOKENS)
    YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_stack_print (yytype_int16 *bottom, yytype_int16 *top)
#else
static void
yy_stack_print (bottom, top)
    yytype_int16 *bottom;
    yytype_int16 *top;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (; bottom <= top; ++bottom)
    YYFPRINTF (stderr, " %d", *bottom);
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)				\
do {								\
  if (yydebug)							\
    yy_stack_print ((Bottom), (Top));				\
} while (YYID (0))


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_reduce_print (YYSTYPE *yyvsp, int yyrule)
#else
static void
yy_reduce_print (yyvsp, yyrule)
    YYSTYPE *yyvsp;
    int yyrule;
#endif
{
  int yynrhs = yyr2[yyrule];
  int yyi;
  unsigned long int yylno = yyrline[yyrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
	     yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      fprintf (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr, yyrhs[yyprhs[yyrule] + yyi],
		       &(yyvsp[(yyi + 1) - (yynrhs)])
		       		       );
      fprintf (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug)				\
    yy_reduce_print (yyvsp, Rule); \
} while (YYID (0))

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef	YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif



#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static YYSIZE_T
yystrlen (const char *yystr)
#else
static YYSIZE_T
yystrlen (yystr)
    const char *yystr;
#endif
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static char *
yystpcpy (char *yydest, const char *yysrc)
#else
static char *
yystpcpy (yydest, yysrc)
    char *yydest;
    const char *yysrc;
#endif
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
	switch (*++yyp)
	  {
	  case '\'':
	  case ',':
	    goto do_not_strip_quotes;

	  case '\\':
	    if (*++yyp != '\\')
	      goto do_not_strip_quotes;
	    /* Fall through.  */
	  default:
	    if (yyres)
	      yyres[yyn] = *yyp;
	    yyn++;
	    break;

	  case '"':
	    if (yyres)
	      yyres[yyn] = '\0';
	    return yyn;
	  }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

/* Copy into YYRESULT an error message about the unexpected token
   YYCHAR while in state YYSTATE.  Return the number of bytes copied,
   including the terminating null byte.  If YYRESULT is null, do not
   copy anything; just return the number of bytes that would be
   copied.  As a special case, return 0 if an ordinary "syntax error"
   message will do.  Return YYSIZE_MAXIMUM if overflow occurs during
   size calculation.  */
static YYSIZE_T
yysyntax_error (char *yyresult, int yystate, int yychar)
{
  int yyn = yypact[yystate];

  if (! (YYPACT_NINF < yyn && yyn <= YYLAST))
    return 0;
  else
    {
      int yytype = YYTRANSLATE (yychar);
      YYSIZE_T yysize0 = yytnamerr (0, yytname[yytype]);
      YYSIZE_T yysize = yysize0;
      YYSIZE_T yysize1;
      int yysize_overflow = 0;
      enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
      char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
      int yyx;

# if 0
      /* This is so xgettext sees the translatable formats that are
	 constructed on the fly.  */
      YY_("syntax error, unexpected %s");
      YY_("syntax error, unexpected %s, expecting %s");
      YY_("syntax error, unexpected %s, expecting %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s");
# endif
      char *yyfmt;
      char const *yyf;
      static char const yyunexpected[] = "syntax error, unexpected %s";
      static char const yyexpecting[] = ", expecting %s";
      static char const yyor[] = " or %s";
      char yyformat[sizeof yyunexpected
		    + sizeof yyexpecting - 1
		    + ((YYERROR_VERBOSE_ARGS_MAXIMUM - 2)
		       * (sizeof yyor - 1))];
      char const *yyprefix = yyexpecting;

      /* Start YYX at -YYN if negative to avoid negative indexes in
	 YYCHECK.  */
      int yyxbegin = yyn < 0 ? -yyn : 0;

      /* Stay within bounds of both yycheck and yytname.  */
      int yychecklim = YYLAST - yyn + 1;
      int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
      int yycount = 1;

      yyarg[0] = yytname[yytype];
      yyfmt = yystpcpy (yyformat, yyunexpected);

      for (yyx = yyxbegin; yyx < yyxend; ++yyx)
	if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
	  {
	    if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
	      {
		yycount = 1;
		yysize = yysize0;
		yyformat[sizeof yyunexpected - 1] = '\0';
		break;
	      }
	    yyarg[yycount++] = yytname[yyx];
	    yysize1 = yysize + yytnamerr (0, yytname[yyx]);
	    yysize_overflow |= (yysize1 < yysize);
	    yysize = yysize1;
	    yyfmt = yystpcpy (yyfmt, yyprefix);
	    yyprefix = yyor;
	  }

      yyf = YY_(yyformat);
      yysize1 = yysize + yystrlen (yyf);
      yysize_overflow |= (yysize1 < yysize);
      yysize = yysize1;

      if (yysize_overflow)
	return YYSIZE_MAXIMUM;

      if (yyresult)
	{
	  /* Avoid sprintf, as that infringes on the user's name space.
	     Don't have undefined behavior even if the translation
	     produced a string with the wrong number of "%s"s.  */
	  char *yyp = yyresult;
	  int yyi = 0;
	  while ((*yyp = *yyf) != '\0')
	    {
	      if (*yyp == '%' && yyf[1] == 's' && yyi < yycount)
		{
		  yyp += yytnamerr (yyp, yyarg[yyi++]);
		  yyf += 2;
		}
	      else
		{
		  yyp++;
		  yyf++;
		}
	    }
	}
      return yysize;
    }
}
#endif /* YYERROR_VERBOSE */


/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep)
#else
static void
yydestruct (yymsg, yytype, yyvaluep)
    const char *yymsg;
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  YYUSE (yyvaluep);

  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  switch (yytype)
    {

      default:
	break;
    }
}


/* Prevent warnings from -Wmissing-prototypes.  */

#ifdef YYPARSE_PARAM
#if defined __STDC__ || defined __cplusplus
int yyparse (void *YYPARSE_PARAM);
#else
int yyparse ();
#endif
#else /* ! YYPARSE_PARAM */
#if defined __STDC__ || defined __cplusplus
int yyparse (void);
#else
int yyparse ();
#endif
#endif /* ! YYPARSE_PARAM */



/* The look-ahead symbol.  */
int yychar;

/* The semantic value of the look-ahead symbol.  */
YYSTYPE yylval;

/* Number of syntax errors so far.  */
int yynerrs;



/*----------.
| yyparse.  |
`----------*/

#ifdef YYPARSE_PARAM
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void *YYPARSE_PARAM)
#else
int
yyparse (YYPARSE_PARAM)
    void *YYPARSE_PARAM;
#endif
#else /* ! YYPARSE_PARAM */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void)
#else
int
yyparse ()

#endif
#endif
{
  
  int yystate;
  int yyn;
  int yyresult;
  /* Number of tokens to shift before error messages enabled.  */
  int yyerrstatus;
  /* Look-ahead token as an internal (translated) token number.  */
  int yytoken = 0;
#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

  /* Three stacks and their tools:
     `yyss': related to states,
     `yyvs': related to semantic values,
     `yyls': related to locations.

     Refer to the stacks thru separate pointers, to allow yyoverflow
     to reallocate them elsewhere.  */

  /* The state stack.  */
  yytype_int16 yyssa[YYINITDEPTH];
  yytype_int16 *yyss = yyssa;
  yytype_int16 *yyssp;

  /* The semantic value stack.  */
  YYSTYPE yyvsa[YYINITDEPTH];
  YYSTYPE *yyvs = yyvsa;
  YYSTYPE *yyvsp;



#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  YYSIZE_T yystacksize = YYINITDEPTH;

  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;


  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY;		/* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */

  yyssp = yyss;
  yyvsp = yyvs;

  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
	/* Give user a chance to reallocate the stack.  Use copies of
	   these so that the &'s don't force the real ones into
	   memory.  */
	YYSTYPE *yyvs1 = yyvs;
	yytype_int16 *yyss1 = yyss;


	/* Each stack pointer address is followed by the size of the
	   data in use in that stack, in bytes.  This used to be a
	   conditional around just the two extra args, but that might
	   be undefined if yyoverflow is a macro.  */
	yyoverflow (YY_("memory exhausted"),
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),

		    &yystacksize);

	yyss = yyss1;
	yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
	goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
	yystacksize = YYMAXDEPTH;

      {
	yytype_int16 *yyss1 = yyss;
	union yyalloc *yyptr =
	  (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
	if (! yyptr)
	  goto yyexhaustedlab;
	YYSTACK_RELOCATE (yyss);
	YYSTACK_RELOCATE (yyvs);

#  undef YYSTACK_RELOCATE
	if (yyss1 != yyssa)
	  YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;


      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
		  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
	YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     look-ahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to look-ahead token.  */
  yyn = yypact[yystate];
  if (yyn == YYPACT_NINF)
    goto yydefault;

  /* Not known => get a look-ahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid look-ahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = YYLEX;
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yyn == 0 || yyn == YYTABLE_NINF)
	goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the look-ahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token unless it is eof.  */
  if (yychar != YYEOF)
    yychar = YYEMPTY;

  yystate = yyn;
  *++yyvsp = yylval;

  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     `$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 2:
#line 116 "kc.y"
    { 
	      SetupTableMan(); 
	    }
    break;

  case 3:
#line 120 "kc.y"
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
	    }
    break;

  case 9:
#line 191 "kc.y"
    { if ((yyvsp[(3) - (3)].strnum).flag==1)
		NewConstant((yyvsp[(1) - (3)].name), (yyvsp[(3) - (3)].strnum).numb);
	      else {
		NewStrConst((yyvsp[(1) - (3)].name), (yyvsp[(3) - (3)].strnum).name);
	      }
	    }
    break;

  case 10:
#line 198 "kc.y"
    { temp=TreeEval((yyvsp[(1) - (2)].tree));
	      if (TreeGetError()==NoEval) 
		fprintf(stderr, "Unable to evaluate expression in line %d.\n", lineno);
	      else {
		(yyval.strnum).flag=1;
		(yyval.strnum).numb=temp;
	      }
	      TreeKill((yyvsp[(1) - (2)].tree));
	    }
    break;

  case 11:
#line 208 "kc.y"
    { (yyval.strnum).flag=2;
	      strcpy((yyval.strnum).name, (yyvsp[(2) - (4)].name));
	    }
    break;

  case 14:
#line 214 "kc.y"
    { NewPrintVar((yyvsp[(1) - (1)].name)); }
    break;

  case 15:
#line 216 "kc.y"
    { NewPrintConc((yyvsp[(1) - (1)].comp).name, (yyvsp[(1) - (1)].comp).charge); }
    break;

  case 17:
#line 219 "kc.y"
    { temp=TreeEval((yyvsp[(3) - (13)].tree));
	      if (TreeGetError()==NoEval)
		fprintf(stderr, "WARNING: Could not evaluate expr in line %d\n", lineno);
	      else {
		NewParameter((yyvsp[(1) - (13)].name), temp);
	      }
	      TreeKill((yyvsp[(3) - (13)].tree)); 
	      temp=TreeEval((yyvsp[(5) - (13)].tree));
	      if (TreeGetError()==NoEval)
		fprintf(stderr, "WARNING: Could not evaluate expr in line %d\n", lineno);
	      else
		NewDeltaParam((yyvsp[(1) - (13)].name), temp);
	      TreeKill((yyvsp[(5) - (13)].tree));
	      temp1=TreeEval((yyvsp[(7) - (13)].tree));
	      if (TreeGetError()==NoEval)
		fprintf(stderr, "WARNING: Could not evaluate expr in line %d\n" , lineno);
	      TreeKill((yyvsp[(7) - (13)].tree));
	      temp2=TreeEval((yyvsp[(11) - (13)].tree));
	      if (TreeGetError()==NoEval)
		fprintf(stderr, "WARNING: Could not evaluate expr in line %d\n" , lineno);
	      TreeKill((yyvsp[(11) - (13)].tree));
	      NewLowHighPrefParam((yyvsp[(1) - (13)].name), temp, temp1, temp2);
	      temp=TreeEval((yyvsp[(13) - (13)].tree));
	      if (TreeGetError()==NoEval)
		fprintf(stderr, "WARNING: Could not evaluate expr in line %d\n", lineno);
	      TreeKill((yyvsp[(13) - (13)].tree));
	      NewDirectForParam((yyvsp[(1) - (13)].name), (int)temp);
	      temp=TreeEval((yyvsp[(9) - (13)].tree));
	      if (TreeGetError()==NoEval)
		fprintf(stderr, "WARNING: Could not evaluate expr in line %d\n", lineno);
	      TreeKill((yyvsp[(9) - (13)].tree));
	      NewDeltaParam((yyvsp[(1) - (13)].name), temp);
	    }
    break;

  case 18:
#line 253 "kc.y"
    { temp=TreeEval((yyvsp[(3) - (13)].tree));
	      if (TreeGetError()==NoEval)
		fprintf(stderr, "WARNING: Could not evaluate expr in line %\n", lineno);
	      else
		NewParamConc((yyvsp[(1) - (13)].comp).name, (yyvsp[(1) - (13)].comp).charge, temp);
	      temp=TreeEval((yyvsp[(5) - (13)].tree));
	      if (TreeGetError()==NoEval)
		fprintf(stderr, "WARNING: Could not evaluate expr in line %\n", lineno);
	      TreeKill((yyvsp[(3) - (13)].tree));
	      TreeKill((yyvsp[(5) - (13)].tree));
	      temp1=TreeEval((yyvsp[(7) - (13)].tree));
	      if (TreeGetError()==NoEval)
		fprintf(stderr, "WARNING: Could not evaluate expr in line %d\n" , lineno);
	      TreeKill((yyvsp[(7) - (13)].tree));
	      temp2=TreeEval((yyvsp[(11) - (13)].tree));
	      if (TreeGetError()==NoEval)
		fprintf(stderr, "WARNING: Could not evaluate expr in line %d\n" , lineno);
	      TreeKill((yyvsp[(11) - (13)].tree));
	      NewLowHighPrefConc((yyvsp[(1) - (13)].comp).name, (yyvsp[(1) - (13)].comp).charge, temp, temp1, temp2);
	      temp=TreeEval((yyvsp[(13) - (13)].tree));
	      if (TreeGetError()==NoEval)
		fprintf(stderr, "WARNING: Could not evaluate expr in line %d\n", lineno);
	      TreeKill((yyvsp[(9) - (13)].tree));
	      NewDirectForConc((yyvsp[(1) - (13)].comp).name, (yyvsp[(1) - (13)].comp).charge, (int)temp);
	      temp=TreeEval((yyvsp[(9) - (13)].tree));
	      if (TreeGetError()==NoEval)
		fprintf(stderr, "WARNING: Could not evaluate expr in line %d\n", lineno);
	      TreeKill((yyvsp[(9) - (13)].tree));
	      NewDeltaConc((yyvsp[(1) - (13)].comp).name, (yyvsp[(1) - (13)].comp).charge, temp);
	    }
    break;

  case 21:
#line 286 "kc.y"
    { NewReaction((int)value); 
	      side=1; 
	    }
    break;

  case 22:
#line 290 "kc.y"
    { 
	      switch ((yyvsp[(5) - (5)].flag)) {
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
	    }
    break;

  case 24:
#line 306 "kc.y"
    { NewExpr((yyvsp[(1) - (5)].name), (yyvsp[(4) - (5)].tree));
	      NewDynVar((yyvsp[(1) - (5)].name));
	      TreeKill((yyvsp[(4) - (5)].tree));
	    }
    break;

  case 29:
#line 315 "kc.y"
    { if (strcmp((yyvsp[(1) - (6)].name), "c")!=0) 
		fprintf(stderr, "c expected!\n");
	      else {
		temp=TreeEval((yyvsp[(6) - (6)].tree));
		if (TreeGetError()==NoEval) 
		  fprintf(stderr, "Unable to evaluate expression in line %d.\n", lineno);
		else
		  NewPowerConst(GetCurrentReaction(), (yyvsp[(3) - (6)].comp).name, (yyvsp[(3) - (6)].comp).charge, temp, side); 
	      }; /* else */
	      TreeKill((yyvsp[(6) - (6)].tree));
	    }
    break;

  case 30:
#line 327 "kc.y"
    { if (strcmp((yyvsp[(1) - (3)].name), "k")==0) {  
		NewRateConst(GetCurrentReaction(), -1, (yyvsp[(3) - (3)].tree)); 
	      } /* if */
	      else 
		if (strcmp((yyvsp[(1) - (3)].name), "v")==0) 
		  NewRateExpr(GetCurrentReaction(), -1, (yyvsp[(3) - (3)].tree));
		else
		  fprintf(stderr, "Syntax error: v or k expected in line %d.\n", lineno);
	      TreeKill((yyvsp[(3) - (3)].tree));               
	    }
    break;

  case 32:
#line 339 "kc.y"
    { if (strcmp((yyvsp[(1) - (3)].name), "K")!=0) 
		fprintf(stderr, "Syntax error: K expected in line %d.\n", lineno);
	      else {
		NewRateConst(GetCurrentReaction(), 0, (yyvsp[(3) - (3)].tree));
	      }; /* else */
	      TreeKill((yyvsp[(3) - (3)].tree)); 
	    }
    break;

  case 33:
#line 346 "kc.y"
    { if (GetReactKind(GetCurrentReaction())==bi) (void) fprintf(stderr, "The reaction in line %d is two-ways.\n", lineno);
	    }
    break;

  case 34:
#line 349 "kc.y"
    { if (GetReactKind(GetCurrentReaction())!=bi)   
		  (void) fprintf(stderr, "The reaction in line %d is a >one-way< reaction or >equilibrium.<\n", lineno);
		else {
		  if (strcmp((yyvsp[(1) - (4)].name), "k")==0)  
		    NewRateConst(GetCurrentReaction(), 1, (yyvsp[(3) - (4)].tree));  
		  else
		    if (strcmp((yyvsp[(1) - (4)].name), "v")==0) 
		      NewRateExpr(GetCurrentReaction(), 1, (yyvsp[(3) - (4)].tree));
		    else
		      fprintf(stderr, "Syntax error: v or k expected in line %d.\n", lineno);
	      }; /* else */
	      TreeKill((yyvsp[(3) - (4)].tree));
	    }
    break;

  case 35:
#line 362 "kc.y"
    { (yyval.flag)='>'; }
    break;

  case 36:
#line 363 "kc.y"
    { (yyval.flag)='<'; }
    break;

  case 37:
#line 364 "kc.y"
    { (yyval.flag)='='; }
    break;

  case 38:
#line 366 "kc.y"
    { SpecieInReaction(GetCurrentReaction(), (yyvsp[(2) - (2)].comp).name, (yyvsp[(2) - (2)].comp).charge);
	      NewCoeff(GetCurrentReaction(), (yyvsp[(2) - (2)].comp).name, (yyvsp[(2) - (2)].comp).charge, (yyvsp[(1) - (2)].dval), side);
	      NewSpecie((yyvsp[(2) - (2)].comp).name, (yyvsp[(2) - (2)].comp).charge);
	    }
    break;

  case 39:
#line 371 "kc.y"
    { SpecieInReaction(GetCurrentReaction(), (yyvsp[(4) - (4)].comp).name, (yyvsp[(4) - (4)].comp).charge);
	      NewCoeff(GetCurrentReaction(), (yyvsp[(4) - (4)].comp).name, (yyvsp[(4) - (4)].comp).charge, (yyvsp[(3) - (4)].dval), side);
	      NewSpecie((yyvsp[(4) - (4)].comp).name, (yyvsp[(4) - (4)].comp).charge);
	    }
    break;

  case 40:
#line 376 "kc.y"
    { 
	      (void) strcpy((yyval.comp).name, (yyvsp[(1) - (2)].name));
	      (yyval.comp).charge=(yyvsp[(2) - (2)].dval);
	      (yyval.comp).concs=no;
	    }
    break;

  case 41:
#line 382 "kc.y"
    { (yyval.dval)=(yyvsp[(2) - (3)].dval); }
    break;

  case 42:
#line 383 "kc.y"
    { 
		(yyval.dval)=0.0;
	      }
    break;

  case 43:
#line 387 "kc.y"
    { charge=FLT_MAX;
	      (yyval.dval)=FLT_MAX;
	    }
    break;

  case 44:
#line 391 "kc.y"
    { if ((yyvsp[(1) - (2)].oper)=='+') {
		(yyval.dval)=(yyvsp[(2) - (2)].dval);
		charge=(yyvsp[(2) - (2)].dval);
	      } /* if */
	      else {
		charge=-(yyvsp[(2) - (2)].dval);
		(yyval.dval)=-(yyvsp[(2) - (2)].dval);
	      }; /* else */
	    }
    break;

  case 45:
#line 401 "kc.y"
    { if ((yyvsp[(1) - (1)].oper)=='+')
		(yyval.dval)=1.0;
	      else
		(yyval.dval)=-1.0;
	    }
    break;

  case 46:
#line 406 "kc.y"
    { (yyval.oper) = '+'; }
    break;

  case 47:
#line 407 "kc.y"
    { (yyval.oper) = '-'; }
    break;

  case 48:
#line 409 "kc.y"
    { coeff=value;
	      (yyval.dval)=(yyvsp[(1) - (1)].dval);  
	    }
    break;

  case 49:
#line 413 "kc.y"
    { coeff=1.0;
	      (yyval.dval)=1.0;
	    }
    break;

  case 52:
#line 419 "kc.y"
    { temp=TreeEval((yyvsp[(4) - (5)].tree)); 
	      switch ((yyvsp[(2) - (5)].flag)) {
	      case 1: 
		if (TreeGetError()==NoEval) 
		  fprintf(stderr, "Unable to evaluate expression.\n");
		else  
		  NewBeginConc((yyvsp[(1) - (5)].comp).name, (yyvsp[(1) - (5)].comp).charge, temp); 
		break;
	      case 0: 
		NewConstraint((yyvsp[(1) - (5)].comp).name, (yyvsp[(1) - (5)].comp).charge, (yyvsp[(4) - (5)].tree));
		break;
	      };
	      TreeKill((yyvsp[(4) - (5)].tree));
	    }
    break;

  case 53:
#line 434 "kc.y"
    { temp=TreeEval((yyvsp[(6) - (7)].tree));
	      NewSpecConst((yyvsp[(3) - (7)].comp).name, (yyvsp[(3) - (7)].comp).charge, (yyvsp[(1) - (7)].name), temp);
	      if (GetError()==NotFound) 
		if ((yyvsp[(3) - (7)].comp).charge==0.0) 
		  NewDynVarConst((yyvsp[(3) - (7)].comp).name, (yyvsp[(1) - (7)].name), temp);
	      TreeKill((yyvsp[(6) - (7)].tree));
	    }
    break;

  case 54:
#line 442 "kc.y"
    { temp=TreeEval((yyvsp[(4) - (5)].tree));
	      if (TreeGetError()==NoEval)
		yyerror("Unable to evaluate expression");
	      else {
		if ((yyvsp[(2) - (5)].flag)==1) {
		  NewInitValue((yyvsp[(1) - (5)].name), temp);
		  TreeKill((yyvsp[(4) - (5)].tree));
		} else
		  yyerror("(0) expected");
	      }
	    }
    break;

  case 55:
#line 454 "kc.y"
    { (yyval.tree)=TreeCreate();
	      TreeCpy((yyval.tree), (yyvsp[(2) - (2)].tree));
	      TreeSign((yyval.tree));
	      TreeKill((yyvsp[(2) - (2)].tree));
	    }
    break;

  case 56:
#line 460 "kc.y"
    { (yyval.tree)=TreeCreate(); 
	      TreeAdd((yyvsp[(1) - (3)].tree), (yyvsp[(3) - (3)].tree));
	      TreeCpy((yyval.tree), (yyvsp[(1) - (3)].tree));
	      TreeKill((yyvsp[(3) - (3)].tree));
	      TreeKill((yyvsp[(1) - (3)].tree));
	    }
    break;

  case 57:
#line 467 "kc.y"
    { (yyval.tree)=TreeCreate(); 
	      TreeSub((yyvsp[(1) - (3)].tree), (yyvsp[(3) - (3)].tree));
	      TreeCpy((yyval.tree), (yyvsp[(1) - (3)].tree));
	      TreeKill((yyvsp[(3) - (3)].tree));
	      TreeKill((yyvsp[(1) - (3)].tree));
	    }
    break;

  case 58:
#line 474 "kc.y"
    {
	      (yyval.tree)=TreeCreate();
	      TreeMul((yyvsp[(1) - (3)].tree), (yyvsp[(3) - (3)].tree));
	      TreeCpy((yyval.tree), (yyvsp[(1) - (3)].tree));
	      TreeKill((yyvsp[(1) - (3)].tree));
	      TreeKill((yyvsp[(3) - (3)].tree));
	    }
    break;

  case 59:
#line 482 "kc.y"
    {
	      (yyval.tree)=TreeCreate();
	      TreeDiv((yyvsp[(1) - (3)].tree), (yyvsp[(3) - (3)].tree));
	      TreeCpy((yyval.tree), (yyvsp[(1) - (3)].tree));
	      TreeKill((yyvsp[(1) - (3)].tree));
	      TreeKill((yyvsp[(3) - (3)].tree));
	    }
    break;

  case 60:
#line 490 "kc.y"
    { (yyval.tree)=TreeCreate();
	      TreeCpy((yyval.tree), (yyvsp[(1) - (3)].tree));
	      TreePow((yyval.tree), (yyvsp[(3) - (3)].tree));
	      TreeKill((yyvsp[(1) - (3)].tree));
	      TreeKill((yyvsp[(3) - (3)].tree));
	    }
    break;

  case 61:
#line 497 "kc.y"
    { (yyval.tree)=TreeCreate();
	      TreeCpy((yyval.tree), (yyvsp[(2) - (3)].tree)); 
	      TreeKill((yyvsp[(2) - (3)].tree));
	    }
    break;

  case 62:
#line 502 "kc.y"
    {
	     (yyval.tree)=TreeCreate();
	     temp=GetConstant((yyvsp[(1) - (1)].name)); 
	     if (GetError()==NotFound) {
	       TreeAssignVar((yyval.tree), (yyvsp[(1) - (1)].name)); 
	       if (strcmp((yyvsp[(1) - (1)].name), "time")==0)
		 NonAutoSystem();
	     } else
	       TreeAssignConst((yyval.tree), temp);
	    }
    break;

  case 63:
#line 513 "kc.y"
    { (yyval.tree)=TreeCreate();
	      TreeAssignConst((yyval.tree), (yyvsp[(1) - (1)].dval));
	    }
    break;

  case 64:
#line 517 "kc.y"
    {  
	      (yyval.tree)=TreeCreate(); 
	      if ((yyvsp[(2) - (2)].flag)==1) {
		temp=GetBeginConc((yyvsp[(1) - (2)].comp).name, (yyvsp[(1) - (2)].comp).charge);
		if (GetError()==NoError)
		  TreeAssignConst((yyval.tree), temp);
		else
		  fprintf(stderr, "[%s(%e)] not found.\n", (yyvsp[(1) - (2)].comp).name, (yyvsp[(1) - (2)].comp).charge);
		flag='1';
	      } /* if */
	      else { 
		flag='0';
		RenameSpec(name, (yyvsp[(1) - (2)].comp).name, (yyvsp[(1) - (2)].comp).charge);
		TreeAssignVar((yyval.tree), name);
	      }; /* else */
	    }
    break;

  case 65:
#line 534 "kc.y"
    { (yyval.tree)=TreeCreate();
	      TreeCpy((yyval.tree), (yyvsp[(3) - (4)].tree));
	      TreeApplyFunc(&(yyval.tree), (yyvsp[(1) - (4)].func));
	      TreeKill((yyvsp[(3) - (4)].tree));
	    }
    break;

  case 66:
#line 539 "kc.y"
    { (yyval.func)=Exp; }
    break;

  case 67:
#line 540 "kc.y"
    { (yyval.func)=Ln; }
    break;

  case 68:
#line 541 "kc.y"
    { (yyval.func)=Log; }
    break;

  case 69:
#line 542 "kc.y"
    { (yyval.func)=Sin; }
    break;

  case 70:
#line 543 "kc.y"
    { (yyval.func)=Cos; }
    break;

  case 71:
#line 544 "kc.y"
    { (yyval.func)=Tan; }
    break;

  case 72:
#line 545 "kc.y"
    { (yyval.func)=Sinh; }
    break;

  case 73:
#line 546 "kc.y"
    { (yyval.func)=Cosh; }
    break;

  case 74:
#line 547 "kc.y"
    { (yyval.func)=Tanh; }
    break;

  case 75:
#line 548 "kc.y"
    { (yyval.func)=Asin; }
    break;

  case 76:
#line 549 "kc.y"
    { (yyval.func)=Acos; }
    break;

  case 77:
#line 550 "kc.y"
    { (yyval.func)=Atan; }
    break;

  case 78:
#line 551 "kc.y"
    { (yyval.func)=Asinh; }
    break;

  case 79:
#line 552 "kc.y"
    { (yyval.func)=Acosh; }
    break;

  case 80:
#line 553 "kc.y"
    { (yyval.func)=Atanh; }
    break;

  case 81:
#line 555 "kc.y"
    { 
	      (void) strcpy((yyval.comp).name, (yyvsp[(2) - (3)].comp).name);
	      (yyval.comp).charge=(yyvsp[(2) - (3)].comp).charge;
	      (yyval.comp).concs=ord;
	    }
    break;

  case 82:
#line 560 "kc.y"
    { (yyval.flag)=0; }
    break;

  case 83:
#line 562 "kc.y"
    { (yyval.flag)=1; }
    break;


/* Line 1267 of yacc.c.  */
#line 2282 "y.tab.c"
      default: break;
    }
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;


  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*------------------------------------.
| yyerrlab -- here on detecting error |
`------------------------------------*/
yyerrlab:
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (YY_("syntax error"));
#else
      {
	YYSIZE_T yysize = yysyntax_error (0, yystate, yychar);
	if (yymsg_alloc < yysize && yymsg_alloc < YYSTACK_ALLOC_MAXIMUM)
	  {
	    YYSIZE_T yyalloc = 2 * yysize;
	    if (! (yysize <= yyalloc && yyalloc <= YYSTACK_ALLOC_MAXIMUM))
	      yyalloc = YYSTACK_ALLOC_MAXIMUM;
	    if (yymsg != yymsgbuf)
	      YYSTACK_FREE (yymsg);
	    yymsg = (char *) YYSTACK_ALLOC (yyalloc);
	    if (yymsg)
	      yymsg_alloc = yyalloc;
	    else
	      {
		yymsg = yymsgbuf;
		yymsg_alloc = sizeof yymsgbuf;
	      }
	  }

	if (0 < yysize && yysize <= yymsg_alloc)
	  {
	    (void) yysyntax_error (yymsg, yystate, yychar);
	    yyerror (yymsg);
	  }
	else
	  {
	    yyerror (YY_("syntax error"));
	    if (yysize != 0)
	      goto yyexhaustedlab;
	  }
      }
#endif
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse look-ahead token after an
	 error, discard it.  */

      if (yychar <= YYEOF)
	{
	  /* Return failure if at end of input.  */
	  if (yychar == YYEOF)
	    YYABORT;
	}
      else
	{
	  yydestruct ("Error: discarding",
		      yytoken, &yylval);
	  yychar = YYEMPTY;
	}
    }

  /* Else will try to reuse look-ahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto yyerrorlab;

  /* Do not reclaim the symbols of the rule which action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;	/* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (yyn != YYPACT_NINF)
	{
	  yyn += YYTERROR;
	  if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
	    {
	      yyn = yytable[yyn];
	      if (0 < yyn)
		break;
	    }
	}

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
	YYABORT;


      yydestruct ("Error: popping",
		  yystos[yystate], yyvsp);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  *++yyvsp = yylval;


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#ifndef yyoverflow
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEOF && yychar != YYEMPTY)
     yydestruct ("Cleanup: discarding lookahead",
		 yytoken, &yylval);
  /* Do not reclaim the symbols of the rule which action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
		  yystos[*yyssp], yyvsp);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  /* Make sure YYID is used.  */
  return YYID (yyresult);
}


#line 563 "kc.y"

#include "lex.c"

