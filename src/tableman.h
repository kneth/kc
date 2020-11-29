/************************************************************************** 
  Symbol Table Manager for kc.

  (C) Copyright 1992-1996 by 
  Kenneth Geisshirt (kneth@fatou.ruc.dk)          Keld Nielsen (kn@kiku.dk)
  Dept. of Life Sciences and Chemistry       Dept. of Theoretical Chemistry 
  Roskilde University                              University of Copenhagen
  P.O. Box 260                                         Universitetsparken 5
  4000 Roskilde                                        2100 Copenhagen East
  Denmark                                                           Denmark

  See kc.tex for details.

  This file contains the definition of constants, types and functions need
  by the (symbol) Table manager.

  Last updated: 2 January 1996 by KG
***************************************************************************/

#ifndef _TABLEMAN_
#define _TABLEMAN_

#include "config.h"
#ifndef _PLATFORM_CONVEX_
#include <float.h>
#endif
#include "symbmath.h" 

#define SymSize       250
#define MaxSpec       150 
#define ReactSize     150
#define ConstrainSize  50 
#define DynSize        50
#define ExprSize       50
#define Conc0Default  0.0
#define MaxSpecConst   10
#define MaxParameter   25 
#define MaxPrint       25
#define StrConstSize   25
#define MaxDynVarConst 10

typedef enum {NoError, TooManyPrn, TooManyConst, TooManySpec, SpecAlready, 
		KonstAlready, NonSpec, TooManyReact, WrongDirect, ReactAlready,
		NotFound, TooManyConstrain, TooManyDynVar, TooManyExpr, 
		ExprAlready, TooManyParam, ParamAlready} TableErrors;

typedef struct SpecK {
  int     used; /* 1=yes, 0=no */ 
  char    name[STRING_LENGTH];
  double  value;
} SpecK;

typedef struct  Spec {
  char       name[STRING_LENGTH];  /* name of species */
  double     charge;          /* charge */
  double     conc0;           /* conc. at t=0 */
  int        is_param;        /* used as parameter */
  SpecK      speck[MaxSpecConst];
} Spec;


typedef struct KKonst {
   char      name[STRING_LENGTH];  /* name of constant */
   double    value;           /* value */
} KKonst;


typedef enum {spec, konst} SymtabTag;

typedef union {
  KKonst   k;
  Spec    s;
}  Data;

typedef struct SymEntry {
  SymtabTag  tag;   
  Data       d;
} SymEntry;

typedef enum {uni, bi, equi} Direc;

typedef struct SpecCoeff {
  char    name[STRING_LENGTH];
  double  charge;
  int     in_use;       /* 1=in use; 0=not in use */
  double  pow_const[2]; /* index=0 <=> reactant, index=1 <=> product */ 
  double  coeff[2];     /* do. */
} SpecCoeff;

typedef struct ReactEntry {
  int        react_no;
  SpecCoeff  species[MaxSpec]; 
  Direc      react_tag;
  int        sel1, sel2; /* selx=1 <=> const, selx=2 <=> rate expr */
  Tree       t1, t2;
} ReactEntry; /* reaction matrix */

typedef struct ConstrainEntry {
  char   name[STRING_LENGTH];
  double charge;
  Tree   expr;
} ConstrainEntry;

typedef struct DynEntry {
  char name[STRING_LENGTH];
  SpecK konst[MaxDynVarConst];
} DynEntry;
  
typedef struct ExprEntry {
  char        name[STRING_LENGTH]; 
  Tree        expr;
  double      init_value;
} ExprEntry;

typedef struct ParamEntry {
  int in_use;         /* 0=no, 1=parameter, 2=concentration */
  char name[STRING_LENGTH];
  double charge; 
  double init_value;  /* initial value */
  double delta;       /* increment */
  double low, high, pref; /* other parameters */
  int    direct;      /* direction for cont. */
} ParamEntry;

typedef struct PrnEntry {
  int tag;            /* 1=name, 2=concentration */
  char name[STRING_LENGTH];
  double charge;
} PrnEntry;

typedef struct StrEntry {
  char name[STRING_LENGTH];
  char value[STRING_LENGTH];
} StrEntry;

SymEntry       symtable[SymSize];     /* symbol table */
int            sym_last;
ReactEntry     reacttable[ReactSize]; /* stoichimetric matrix */
int            react_last;
int            current;               /* pointer current reaction */
ConstrainEntry contable[ConstrainSize]; /* constrain table */
int            con_last;
DynEntry       dyntable[DynSize];     /* dynamical variables */
int            dyn_last;
TableErrors    error;                
int            search_i, search_j, search_no;   /* search vars */
ExprEntry      exprtable[ExprSize];   /* expression table */
int            expr_last;
ParamEntry     paramtable[MaxParameter];
int            param_last;
PrnEntry       prntable[MaxPrint];
int            prn_last;
StrEntry       strtable[StrConstSize];
int            str_last;
int            autonom_system; /* 0 if autonomous ODEs, 1 otherwise */

extern void         SetupTableMan(void);
extern TableErrors  GetError(void);
extern void         NewSpecie(char *, double);
extern void         NewConstant(char *, double);
extern void         NewRateConst(int, int, Tree);
extern void         NewRateExpr(int, int, Tree);
extern void         NewBeginConc(char *, double, double);
extern void         SpecieInReaction(int, char *, double);
extern void         NewReaction(int);
extern int          GetCurrentReaction(void);
extern void         AddReactionKind(int, Direc);
extern int          GetReactionNo(int);
extern void         RenameSpec(char *, char *, double);
extern int          NoOfSpec(void);
extern int          GetFirstSpecA(int, char *, double *, double *, int);
extern int          GetNextSpecA(char *, double *, double *, int);
extern int          IsSpecInReact(int, char *, double, double *);
extern int          GetFirstSpecB(char *, double *);
extern int          GetNextSpecB(char *, double *);
extern double       GetCoeffInReact(int, char *, double, int);
extern double       GetPowConstInReact(int, char *, double, int);
extern void         GetSpecNo(int, char *, double *);
extern Direc        GetReactKind(int);
extern void         GetRateConst(int, Direc, int, Tree);
extern void         GetRateExpr(int, Direc, int, Tree);
extern int          NoOfReact(void);
extern void         NewCoeff(int, char *, double, double, int);
extern double       SumCoeff(int, int);
extern double       GetConstant(char *);
extern double       GetBeginConc(char *, double);
extern int          GetSpecNumber(char *, double);
extern void         NewDynVar(char *);
extern int          NumOfDynVar(void);
extern void         GetDynVarNo(int, char *);
extern void         NewExpr(char *, Tree);
extern int          NumOfExpr(void);
extern void         GetExprNo(int, char *, Tree);
extern int          IsVarParameter(char *);
extern void         NewPowerConst(int, char *, double, double, int);
extern void         NewConstraint(char *, double, Tree);
extern int          NumOfConstraint(void);
extern void         GetConstraintNo(int, char *, double *, Tree);
extern int          IsSpecInConstraint(char *, double);
extern void         NewSpecConst(char *, double, char *, double);
extern double       GetSpecConst(char *, double, char *);
extern void         NewParamter(char *, double);
extern void         GetParamNo(int, char *, double *, int *);
extern int          NumOfParameter(void);
extern void         NewParamConc(char *, double, double);
extern int          IsSpecParam(char *, double);
extern void         GetConstantNo(int i, char *name);
extern int          NumOfConstants(void);
extern void         NewDeltaParam(char *, double);
extern void         NewDeltaConc(char *, double, double);
extern void         GetDeltaParam(char *, double *);
extern void         GetDeltaConc(char *, double, double *);
extern void         GetInitParam(char *, double *);
extern void         GetInitConc(char *, double, double *); 
extern void         NewPrintVar(char *);
extern void         NewPrintConc(char *, double);
extern int          NumOfPrint(void);
extern void         GetPrintNo(int, char *, double *, int *);
extern int          IsSpecInPrnList(char *, double, int);
extern void         NewLowHighPrefParam(char *, double, double, double);
extern void         GetLowHighPrefParam(char *, double *, double *, double *);
extern void         NewLowHighPrefConc(char *, double, double, double, double);
extern void         GetLowHighPrefConc(char *, double, double *, double *, 
				       double *);
extern void         NewInitValue(char *, double);
extern double       GetInitValue(char *);
extern int          NoSpecInReacs(char *, double);
extern void         NewDirectForParam(char *, int);
extern int          GetDirectForParam(char *);
extern void         NewDirectForConc(char *, double, int);
extern int          GetDirectForConc(char *, double);
extern void         NewStrConst(char *, char *);
extern void         GetStrConst(char *, char *);
extern int          NumOfStrConst(void);
extern void         GetStrConstNo(int, char *, char *);
extern void         NonAutoSystem(void);
extern int          IsNonAutoSystem(void);
extern void         NewDynVarConst(char *, char *, double);
extern double       GetDynVarConst(char *, char *);
extern void         GetStocMatrix(int, int, double **);
#endif
