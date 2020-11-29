/************************************************************************** 
  Implementation of the Symbol Table Manager for kc.

  (C) Copyright 1992-1996 by 
  Kenneth Geisshirt (kneth@fatou.ruc.dk)      Keld Nielsen (kn@kin.kiku.dk)
  Dept. of Life Sciences and Chemistry       Dept. of Theoretical Chemistry 
  Roskilde University                              University of Copenhagen
  P.O. Box 260                                         Universitetsparken 5
  4000 Roskilde                                        2100 Copenhagen East
  Denmark                                                           Denmark

  Last updated: 27 March 1996 by KG
***************************************************************************/

#include "config.h"

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "tableman.h"

void SetupTableMan(void) {

  sym_last=0;
  react_last=0;
  con_last=0;
  current=0;
  dyn_last=0;
  expr_last=0;
  param_last=0;
  prn_last=0;
  str_last=0;
  error=NoError;
  autonom_system=0;
}  /* SetupTableMan */

TableErrors GetError(void) {

  return error;
} /* GetError */

int GetCurrentReaction(void) {

  return current;
} /* GetCurrentReaction */

int NoOfReact(void) {

  return react_last;
} /* NoOfReact */

void NewSpecie(char *name, double charge) {

  int i, finished, j;

  if (sym_last==SymSize)  
    error=TooManySpec;
  else {
    finished=0;
    i=0;
    while ((i<sym_last) && (finished==0)) {
      if (symtable[i].tag==spec) {
        if ((strcmp(name, symtable[i].d.s.name)==0) && 
	    (symtable[i].d.s.charge==charge)) {
	  finished=1;
          error=SpecAlready;
        } /* if */
      } /* if */
      i++;
    } /* while */
    if (finished==0) {
      symtable[sym_last].tag=spec;
      (void) strcpy(symtable[sym_last].d.s.name, name);
      symtable[sym_last].d.s.charge=charge;
      symtable[sym_last].d.s.conc0=Conc0Default;
      for(j=0; j<MaxSpecConst; j++)
	symtable[i].d.s.speck[j].used=0;
      symtable[sym_last].d.s.is_param=0;
      sym_last++;
      error=NoError;
    } /* if */
  } /* else */
} /* NewSpecie */

void NewSpecConst(char *name1, double charge, char *name2, double value) {

  int i, j, finished;

  finished=0;
  i=0;
  while ((i<sym_last) && (finished==0)) {
    if (symtable[i].tag==spec) {
      if ((strcmp(symtable[i].d.s.name, name1)==0) && 
	  (symtable[i].d.s.charge==charge)) {
	j=0;
	while ((j<MaxSpecConst) && (finished==0)) {
	  if (symtable[i].d.s.speck[j].used==0) {
	    strcpy(symtable[i].d.s.speck[j].name, name2);
	    symtable[i].d.s.speck[j].value=value;
	    symtable[i].d.s.speck[j].used=1;
	    finished=1;
          } /* if */
	  j++;
        } /* while */
      } /* if */
    } /* if */
    i++;
  } /* while */
  if (finished==1)
    error=NoError;
  else
    error=NotFound;
} /* NewSpecConst */

double GetSpecConst(char *name1, double charge, char *name2) {

  int i, j, finished;
  double temp;

  finished=0;
  i=0;
  while ((i<sym_last) && (finished==0)) {
    if (symtable[i].tag==spec) {
      if ((strcmp(symtable[i].d.s.name, name1)==0) && 
	  (symtable[i].d.s.charge==charge)) {
	j=0;
	while ((j<MaxSpecConst) && (finished==0)) {
	  if (symtable[i].d.s.speck[j].used==1) {
	    if (strcmp(symtable[i].d.s.speck[j].name, name2)==0) {
	      temp=symtable[i].d.s.speck[j].value;
	      finished=1;
            } /* if */
          } /* if */
	  j++;
        } /* while */
      } /* if */
    } /* if */
    i++;
  } /* while */
  if (finished==1)
    error=NoError;
  else
    error=NotFound;
  return temp;
} /* GetSpecConst */

void NewDynVarConst(char *dynvar, char *name, double value) {

  int i, j, finished;

  finished=0;
  i=0;
  while ((i<dyn_last) && (finished==0)) {
    if (strcmp(dyntable[i].name, dynvar)==0)
      j=0;
    while ((j<MaxDynVarConst) && (finished==0)) {
      if (dyntable[i].konst[j].used==0) {
	strcpy(dyntable[i].konst[j].name, name);
	dyntable[i].konst[j].value=value;
	dyntable[i].konst[j].used=1;
	finished=1;
      } /* if */
      j++;
    } /* while */
    i++;
  } /* while */
  if (finished==1)
    error=NoError;
  else
    error=NotFound;
} /* NewDynVarConst */

double GetDynVarConst(char *name1, char *name2) {

  int i, j, finished;
  double temp;

  finished=0;
  i=0;
  while ((i<dyn_last) && (finished==0)) {
    if (strcmp(dyntable[i].name, name1)==0)
      j=0;
    while ((j<MaxDynVarConst) && (finished==0)) {
      if (dyntable[i].konst[j].used==1) {
	if (strcmp(dyntable[i].konst[j].name, name2)==0) {
	  temp=dyntable[i].konst[j].value;
	  finished=1;
	} /* if */
      } /* if */
      j++;
    } /* while */
    i++;
  } /* while */
  if (finished==1)
    error=NoError;
  else
    error=NotFound;
  return temp;
} /* GetDynVarConst */

void NewConstant(char *name, double value) {

  int i, finished;

  error=NoError;
  if (sym_last==SymSize)
    error=TooManyConst;
  else {
    finished=0;
    i=0;
    while ((i<sym_last) && (finished==0)) {
      if (symtable[i].tag==konst) {
	if (strcmp(name, symtable[i].d.k.name)==0) { 
	  finished=1;
          error=KonstAlready;
        } /* if */
      } /* if */
      i++;
    } /* while */
    if (error==NoError) {
      symtable[sym_last].tag=konst;
      (void) strcpy(symtable[sym_last].d.k.name, name);
      symtable[sym_last].d.k.value=value;
      sym_last++;
      error=NoError;
    } /* if */
  } /* else */
} /* NewConstant */

void NewCoeff(int react_no, char *name, double charge, 
	      double coeff, int side) {

  int i, j, finished;

  i=0;
  finished=0;
  while ((i<react_last) && (finished==0)) {
    if (reacttable[i].react_no==react_no) {
      j=0;
      finished=0;  /* not needed!! */
      while ((j<MaxSpec) && (finished==0)) {
	if (reacttable[i].species[j].in_use==1) {
	  if ((reacttable[i].species[j].charge==charge) && 
	      (strcmp(reacttable[i].species[j].name, name)==0)) {
            if (side==1) {
	      reacttable[i].species[j].coeff[0]+=coeff;
	      reacttable[i].species[j].pow_const[0]+=coeff;
            } else {
              reacttable[i].species[j].coeff[1]+=coeff;
              reacttable[i].species[j].pow_const[1]+=coeff;
            } /* if ... else */
	    finished=1;
          } /* if */
          else
	    j++;
        } /* if */
      } /* while */
    } /* if */
    else 
      i++;
  } /* while */
} /* NewCoeff */

double SumCoeff(int react_no, int side) {

  int i, j, finished;
  double res;

  i=0;
  finished=0;
  while ((i<react_last) && (finished==0)) {
    if (reacttable[i].react_no==react_no) {
      j=0;
      res=0.0;
      while (j<MaxSpec) {
	if (reacttable[i].species[j].in_use==1) {
	  if (side==1) 
	    res+=reacttable[i].species[j].coeff[0];
	  else
	    res+=reacttable[i].species[j].coeff[1];
	}
	j++;
      } /* while */
      finished=1;
    } else
      i++;
  } /* while */
  return res;
} /* SumCoeff */

void NewPowerConst(int react_no, char *name, double charge, 
		   double value, int side) {

  int i, j, finished;

  i=0;
  finished=0;
  while ((i<react_last) && (finished==0)) {
    if (reacttable[i].react_no==react_no) {
      j=0;
      finished=0;  /* not needed!! */
      while ((j<MaxSpec) && (finished==0)) {
	if (reacttable[i].species[j].in_use==1) {
	  if ((reacttable[i].species[j].charge==charge) && 
	      (strcmp(reacttable[i].species[j].name, name)==0)) { 
            if (side==1)
	      reacttable[i].species[j].pow_const[0]=value;
            else
              reacttable[i].species[j].pow_const[1]=value;
	    finished=1;
          } /* if */
        } /* if */
	j++;
      } /* while */
    } /* if */
    i++;
  } /* while */
} /* NewPowerConst */

void NewRateConst(int react, int direct, Tree value) {
 
/* For direct: -1 = >, 0 = =, 1 = < */

  int i, finished;

  i=0;
  finished=0;
  while ((i<react_last) && (finished==0)) {
    if (reacttable[i].react_no==react) {
      switch (reacttable[i].react_tag) {
      case uni: 
	if (direct==-1) {
	  reacttable[i].t1=TreeCreate();
	  TreeCpy(reacttable[i].t1, value);
	  reacttable[i].sel1=1;
	  error=NoError;
	} /* if */
	else 
	  error=WrongDirect;
	break;
      case bi: 
	if (direct==-1) {
	  reacttable[i].t1=TreeCreate();
	  TreeCpy(reacttable[i].t1, value);
	  reacttable[i].sel1=1;
	  error=NoError;
	} /* if */
	else if (direct==1) {
          reacttable[i].t2=TreeCreate();
	  TreeCpy(reacttable[i].t2, value);
	  reacttable[i].sel2=1;
	  error=NoError;
	} /* else-if */
	else 
	  error=WrongDirect;
	break;
      case equi: 
	if (direct==0) {
          reacttable[i].t1=TreeCreate();
	  TreeCpy(reacttable[i].t1, value);
	  reacttable[i].sel1=1;	
	  error=NoError;
	} /* if */
	else
	  error=WrongDirect;
	break;
      } /* switch */
      finished=1;
    } /* if */
    i++;
  } /* while */
} /* NewRateConst */

void NewRateExpr(int react, int direct, Tree value) {
 
/* For direct: -1 = >, 0 = =, 1 = < */

  int i, finished;

  i=0;
  finished=0;
  while ((i<react_last) && (finished==0)) {
    if (reacttable[i].react_no==react) {
      switch (reacttable[i].react_tag) {
      case uni: 
	if (direct==-1) {
	  reacttable[i].t1=TreeCreate();
	  TreeCpy(reacttable[i].t1, value);
	  reacttable[i].sel1=2;
	  error=NoError;
	} /* if */
	else 
	  error=WrongDirect;
	break;
      case bi: 
	if (direct==-1) {
	  reacttable[i].t1=TreeCreate();	
	  TreeCpy(reacttable[i].t1, value);
	  reacttable[i].sel1=2;
	  error=NoError;
	} /* if */
	else if (direct==1) {
	  reacttable[i].t2=TreeCreate();
	  TreeCpy(reacttable[i].t2, value);
	  reacttable[i].sel2=2;
	  error=NoError;
	} /* else-if */
	else 
	  error=WrongDirect;
	break;
      case equi :
	error=WrongDirect;
	break;
      } /* switch */
      finished=1;
    } /* if */
    i++;
  } /* while */
} /* NewRateExpr */

void NewBeginConc(char *name, double charge, double value) {

  int i, finished;

  finished=0;
  i=0;
  while ((i<sym_last) && (finished==0)) {
    if (symtable[i].tag==spec) {
      if ((strcmp(name, symtable[i].d.s.name)==0) && 
	  (symtable[i].d.s.charge==charge)) {
        finished=1;
      } /* if */
    } /* if */
    i++;
  } /* while */
  if (finished==1) {
    symtable[i-1].d.s.conc0=value;
    error=NoError;
  } /* if */
  else 
    error=NonSpec;
} /* NewBeginConc */


void NewReaction(int react) {

  int i, j, finished;
  
  if (react_last==ReactSize)
    error=TooManyReact;
  else {
    finished=0;
    i=0;
    while ((i<react_last) && (finished==0)) {
      if (reacttable[i].react_no==react) {
	error=ReactAlready;
	finished=1;
      } /* if */
      else
	i++;
     } /* while */
     if ((i==react_last) && (finished==0)) {
       reacttable[i].react_no=react;
       for(j=0; j<MaxSpec; j++) {
	 reacttable[i].species[j].in_use=0;
         reacttable[i].species[j].coeff[0]=0.0;   /* reset coeff */
	 reacttable[i].species[j].coeff[1]=0.0;
	 reacttable[i].species[j].pow_const[0]=0.0;
	 reacttable[i].species[j].pow_const[1]=0.0;
       } /* for */
       react_last++;
     } /* if */
  } /* else */
  current=react;
} /* NewReaction */

void AddReactionKind(int react, Direc direct) {

  int  finished, i;

  error=NotFound;
  finished=0;
  i=0;
  while ((i<=react_last) && (finished==0)) {
    if (reacttable[i].react_no==react) {
      finished=1;
      reacttable[i].react_tag=direct;
      error=NoError;
    } /* if */
    i++;
  } /* while */
} /* AddReactionKind */


void SpecieInReaction(int react, char *name, double charge) {

  int i, j, finished;

  i=0;
  finished=0;
  while ((i<react_last) && (finished==0)) {
    if (reacttable[i].react_no==react) {
      j=0;
      finished=0;  /* not needed!! */
      while ((j<MaxSpec) && (finished==0)) {
	if (reacttable[i].species[j].in_use==1) {
	  if ((strcmp(reacttable[i].species[j].name, name)==0) && 
	      (reacttable[i].species[j].charge==charge)) {
            finished=1;
            error=SpecAlready;
          } /* if */
	  else
	    j++;
	} /* if */
	else {
	  reacttable[i].species[j].in_use=1;
	  (void) strcpy(reacttable[i].species[j].name, name);
	  reacttable[i].species[j].charge=charge;
	  finished=1;
	} /* else */
      } /* while */
    } /* if */
    else 
      i++;
  } /* while */
} /* SpecieInReaction */

int GetReactNo(int counter) {
 
/****
  Please note:
  No error check is done; the user is assumed to know what he is doing 
****/

 return reacttable[counter].react_no;
} /* GetReactNo */


void RenameSpec(char *rename, char *name, double charge) {

  int i;

/* charge = FLT_MAX => specie is radical */
 
  (void) strcpy(rename, name);
  if (charge==FLT_MAX)
    (void) strcat(rename, "_rad\0");
  else if (charge==0.0)
    /* nothing */ ;
    else {
      (void) strcat(rename, "_\0");
      if (charge>=0.0) 
       for(i=1; i<=(int)charge; i++) 
	 (void) strcat(rename, "p\0");
      else
	for(i=1; i<=-(int)charge; i++)  
	  (void) strcat(rename, "m\0");
    } /* else */
} /* RenameSpec */

int NoOfSpec(void) {

  int i, count;

  count=0;
  for(i=0;i<sym_last;i++) {
    if (symtable[i].tag==spec) 
      count++;
  } /* for */
  return count;
} /* NoOfSpec */

int GetFirstSpecA(int no, char *name, double *charge, 
		  double *coeff, int side) {

  int finished, i, j;

  i=0;
  finished=0;
  while ((i<react_last) && (finished==0)) {
    if (reacttable[i].react_no==no) {
      j=0;
      while ((j<MaxSpec) && (finished==0)) {
	if (reacttable[i].species[j].in_use==1) {
	    finished=1;
	    (void) strcpy(name, reacttable[i].species[j].name);
	    *charge=reacttable[i].species[j].charge;
	    *coeff=reacttable[i].species[j].coeff[side];
        } /* if */
	j++;   
      } /* while */
    } /* if */
    i++;
  } /* while */
  search_i=i-1;
  search_j=j;
  search_no=no;
  return finished;
} /* GetFirstSpecA */ 

int GetNextSpecA(char *name, double *charge, double *coeff, int side) {

  int i, j, finished;

  finished=0;
  i=search_i;
  j=search_j;
      while ((j<MaxSpec) && (finished==0)) {
	if (reacttable[i].species[j].in_use==1) {
	  finished=1;
	  (void) strcpy(name, reacttable[i].species[j].name);
	  *charge=reacttable[i].species[j].charge;
	  *coeff=reacttable[i].species[j].coeff[side];
	} /* if */
	j++;
      } /* while */
  search_j=j;
  return finished;
} /* GetNextSpecA */ 

int IsSpecInReact(int react_no, char *name, double charge, double *coeff) {
  
  int i, j, finished;

  i=0;
  finished=0;
  while ((i<react_last) && (finished==0)) {
    if (reacttable[i].react_no==react_no) {
      j=0;
      while (j<MaxSpec) {
	if (reacttable[i].species[j].in_use==1) {
	  if ((strcmp(name, reacttable[i].species[j].name)==0) && 
	      (reacttable[i].species[j].charge==charge)) {  
	    *coeff=reacttable[i].species[j].coeff[0]
	      -reacttable[i].species[j].coeff[1];
	    finished=1;
          } /* if */
        } /* if */ 
	j++;
      } /* while */
    } /* if */
    i++;
  } /* while */
  return finished;
} /* IsSpecInReact */

double GetCoeffInReact(int react_no, char *name, double charge, int side) {
  
  int i, j;

  i=0;
  while (i<react_last) {
    if (reacttable[i].react_no==react_no) {
      j=0;
      while (j<=MaxSpec) {
	if (reacttable[i].species[j].in_use==1) {
	  if ((strcmp(name, reacttable[i].species[j].name)==0) && 
	      (reacttable[i].species[j].charge==charge)) {
            if (side==1)  
	      return reacttable[i].species[j].coeff[0];
            else
              return reacttable[i].species[j].coeff[1];
          } /* if */
        } /* if */ 
	else
	  j++;
      } /* while */
    } /* if */
    i++;
  } /* while */
  return 0.0;
} /* GetCoeffInReact */

double GetPowConstInReact(int react_no, char *name, double charge, int side) {
  
  int i, j, finished;
  double temp;

  i=0;
  finished=0;
  while ((i<react_last) && (finished==0)) {
    if (reacttable[i].react_no==react_no) {
      j=0;
      while ((j<MaxSpec) && (finished==0)) {
	if (reacttable[i].species[j].in_use==1) {
	  if ((strcmp(name, reacttable[i].species[j].name)==0) && 
	      (reacttable[i].species[j].charge==charge)) {
	    temp=reacttable[i].species[j].pow_const[side];
	    finished=1;
          } /* if */
        } /* if */ 
	j++;
      } /* while */
    } /* if */
    i++;
  } /* while */
  return temp;
} /* GetPowConstInReact */


int GetFirstSpecB(char *name, double *charge) {

  int finished, i;

  i=0;
  finished=0;
  while ((i<sym_last) && (finished==0)) {
    if (symtable[i].tag==spec) {
      finished=1;
      search_i=i;
      (void) strcpy(name, symtable[i].d.s.name);
      *charge=symtable[i].d.s.charge;
    } /* if */
    i++;
  } /* while */
  return finished;
} /* GetFirstSpecB */

int GetNextSpecB(char *name, double *charge) {

  int finished, i;

  i=search_i+1;
  finished=0;
  while ((i<sym_last) && (finished==0)) {
    if (symtable[i].tag==spec) {
      finished=1;
      search_i=i;
      (void) strcpy(name, symtable[i].d.s.name);
      *charge=symtable[i].d.s.charge;
    } /* if */
    i++;
  } /* while */
  return finished;
} /* GetNextSpecB */

void GetSpecNo(int count, char *name, double *charge) {
 
  int i;

  i=0;
  while ((count>0) && (i<sym_last)) {
    if (symtable[i].tag==spec) 
      count--;
    i++;
  } /* while */
  (void) strcpy(name, symtable[i-1].d.s.name);
  *charge=symtable[i-1].d.s.charge;
} /* GetSpecNo */


Direc GetReactKind(int react_no) {

  int i, finished;

  finished=0;
  i=0;
  while ((i<react_last) && (finished==0)) {
    if (reacttable[i].react_no==react_no)  
      finished=1;
    else
      i++;
  } /* while */
  return reacttable[i].react_tag;
} /* GetReactKind */

void GetRateConst(int react_no, Direc direct, int way, Tree res) {

/****
  way is only used when direct==bi, 
  way==1 <=> left-to-right
  way==2 <=> right-to-left
****/

  int i, finished;

  finished=0;
  i=0;
  while ((i<react_last) && (finished==0)) {
    if ((reacttable[i].react_no==react_no) && 
	(reacttable[i].react_tag==direct))  
      finished=1;
    else
      i++;
  } /* while */ 
  if (direct==bi) 
    if (way==1) 
      TreeCpy(res, reacttable[i].t1);
    else
      TreeCpy(res, reacttable[i].t2);
  else
    TreeCpy(res, reacttable[i].t1); 
} /* GetRateConst */     

void GetRateExpr(int react_no, Direc direct, int way, Tree t) {

/****
  way is only used when direct==bi, 
  way==1 <=> left-to-right
  way==2 <=> right-to-left
****/

  int i, finished;

  finished=0;
  i=0;
  while ((i<react_last) && (finished==0)) {
    if ((reacttable[i].react_no==react_no) && 
	(reacttable[i].react_tag==direct))  
      finished=1;
    else
      i++;
  } /* while */ 
  if (direct==bi) 
    if (way==1) 
      TreeCpy(t, reacttable[i].t1);
    else
      TreeCpy(t, reacttable[i].t2);
  else
    TreeCpy(t, reacttable[i].t1); 
} /* GetRateExpr */     

int GetRateKind(int react_no, Direc direct, int way) {

/***** 2=expression, 1=constant *****/

  int i, finished;

  finished=0;
  i=0;
  while ((i<react_last) && (finished==0)) {
    if ((reacttable[i].react_no==react_no) && 
	(reacttable[i].react_tag==direct))
      finished=1;
    else
      i++;
  } /* while */
  if (direct==bi) 
    if (way==1)
      return reacttable[i].sel1;
    else 
      return reacttable[i].sel2;
  else
    return reacttable[i].sel1;
} /* GetRateKind */

double GetConstant(char *name) {

  int i, finished;
  double temp=0.0;

  finished=0;
  i=0;
  while ((i<sym_last) && (finished==0)) {
    if (symtable[i].tag==konst) { 
      if (strcmp(symtable[i].d.k.name, name)==0) {
	finished=1;
	error=NoError;
	temp=symtable[i].d.k.value;
      } /* if */
    } /* if */
    i++;
  } /* while */
  if (finished==0)
    error=NotFound;
  return temp;
} /* GetConstant */
      
double GetBeginConc(char *name, double charge) {

  int i, finished;
  double temp;

  finished=0;
  i=0;
  while ((i<sym_last) && (finished==0)) {
    if (symtable[i].tag==spec) { 
      if ((strcmp(symtable[i].d.s.name, name)==0) && 
	  (symtable[i].d.s.charge==charge)) {
	finished=1;
	error=NoError;
	temp=symtable[i].d.s.conc0;
      } /* if */
    } /* if */
    i++;
  } /* while */
  if (finished==0)
    error=NotFound;
  return temp;
} /* GetBeginConc */

int GetSpecNumber(char *name, double charge) {

  int counter;
  int finished;
  int i;

  counter=0;
  finished=0;
  i=0;
  while ((finished==0) && (i<sym_last)) {
    if (symtable[i].tag==spec) {
      counter++;
      if ((strcmp(symtable[i].d.s.name, name)==0) && 
	  (symtable[i].d.s.charge==charge))
	finished=1;
    } /* if  */
    i++;
  } /* while */
  return counter;
} /* GetSpecNumber */

void NewDynVar(char *name) {

  if (dyn_last==DynSize)  
    error=TooManyDynVar;
  else { 
    strcpy(dyntable[dyn_last].name, name);
    dyn_last++;
    error=NoError;
  } /* else */
} /* NewDynVar */

int NumOfDynVar(void) {

  return dyn_last;
} /* NumOfDynVar */

void GetDynVarNo(int i, char *name) {

  if ((i-1)>=dyn_last) 
    error=NotFound;
  else {
    strcpy(name, dyntable[i-1].name);
    error=NoError;
  } /* else */
} /* GetDynVarNo */

void NewExpr(char *name, Tree t) {

  int i, found;

  if (expr_last>=ExprSize) 
    error=TooManyExpr;
  else {
    i=0;
    found=0;
    while ((i<expr_last) && (found==0)) {
      if (strcmp(exprtable[i].name, name)==0) {
	error=ExprAlready;
	found=1;
      } /* if */
      else 
	i++;
    } /* while */
    if (found==0) {
      strcpy(exprtable[expr_last].name, name);
      exprtable[expr_last].expr=TreeCreate();
      TreeCpy(exprtable[expr_last].expr, t);
      expr_last++;
      error=NoError;
    } /* if */
  } /* else */
} /* NewExpr */
    
int NumOfExpr(void) {

  return expr_last;
} /* NumOfExpr */

void GetExprNo(int no, char *name, Tree t) {

  if (no>expr_last) {
    error=NotFound;
    t=NULL;
  } /* if */
  else {
    error=NoError;
    strcpy(name, exprtable[no-1].name);
    TreeCpy(t, exprtable[no-1].expr);
  } /* else */
} /* GetExprNo */

int IsVarParameter(char *name) {

/****
  1 : Variable is parameter   
  0 : otherwise
****/

  int found, i;

  i=0;
  found=0;
  while ((i<=expr_last) && (found==0)) {
    if (strcmp(name, exprtable[i].name)==0)
      found=1;
    else
      i++;
  } /* while */
  return found;
} /* IsVarParameter */

void NewConstraint(char *name, double charge, Tree t) {

  if (con_last==ConstrainSize)  
    error=TooManyConstrain;
  else {
    contable[con_last].expr=TreeCreate();
    strcpy(contable[con_last].name, name);
    TreeCpy(contable[con_last].expr, t);
    contable[con_last].charge=charge;
    con_last++;
    error=NoError;
  } /* else */
} /* NewConstaint */

int NumOfConstraint(void) {

  return con_last;
} /* NumOfConstraint */

void GetConstraintNo(int no, char *name, double *charge, Tree t) {

  if (no>con_last) { 
    error=NotFound;
    t=NULL;
  } /* if */
  else {
    strcpy(name, contable[no-1].name);
    *charge=contable[no-1].charge;
    error=NoError;
    TreeCpy(t, contable[no-1].expr);
  } /* else */
} /* GetConstraintNo */
   
int IsSpecInConstraint(char *name, double charge) {

/****
   return value : 0  - not found
                  >0 - constraint number
****/

  int i, finished, found;

  found=0;
  finished=0;
  i=0;

  while ((i<con_last) && (finished==0)) {
    if ((strcmp(name, contable[i].name)==0) && (charge==contable[i].charge)) {
      finished=1;
      found=1;
    } /* if */
    i++;
  } /* while */
  if (found==1) 
    return i;
  else
    return 0;
} /* IsSpecInConstraint */

void NewParameter(char *name, double init) {

  int i, finished;

  if (param_last==(MaxParameter-1)) 
    error=TooManyParam;
  else {
    i=0;
    finished=0;
    while ((i<param_last) && (finished==0)) {
      if (paramtable[i].in_use==1) { 
	if (strcmp(paramtable[i].name, name)==0) {
	  error=ParamAlready;
	  finished=1;
	} /* if */
      }
      i++;
    } /* while */
    if ((i==param_last) && (finished==0)) {
      strcpy(paramtable[i].name, name);
      paramtable[i].init_value=init;
      paramtable[i].in_use=1;
      param_last++;
      error=NoError;
    }
  } /* else */
} /* NewParameter */

void NewDeltaParam(char *name, double delta) {

  int i, finished;

  finished=0;
  i=0;
  while ((i<param_last) && (finished==0)) {
    if (paramtable[i].in_use==1) {
      if (strcmp(paramtable[i].name, name)==0) { 
	paramtable[i].delta=delta;
	finished=1;
      }
    } /* if */
    i++;
  } /* while */
} /* NewDeltaParam */

void NewDeltaConc(char *name, double charge, double delta) {

  int i, finished;
  
  finished=0;
  i=0;
  while ((i<param_last) && (finished==0)) {
    if (paramtable[i].in_use==2) {
      if ((strcmp(paramtable[i].name, name)==0) && 
	  (paramtable[i].charge==charge)) {
	paramtable[i].delta=delta;
	finished=1;
      } /* if */
    } /* if */
    i++;
  } /* while */
} /* NewDeltaConc */

void NewParamConc(char *name, double charge, double init) {

  int i, j, finished;

  if (param_last==(MaxParameter-1))
    error=TooManyParam;
  else {
    i=0;
    finished=0;
    while ((i<param_last) && (finished==0)) {
      if (paramtable[i].in_use==2) {
        if ((strcmp(paramtable[i].name, name)==0) && 
	    paramtable[i].charge==charge)
          error=ParamAlready;
      } 
      i++;
    } /* while */
    if ((i==param_last) && (finished==0)) {
      strcpy(paramtable[i].name, name);
      paramtable[i].charge=charge;
      paramtable[i].init_value=init;
      paramtable[i].in_use=2;
      param_last++;
      i=0;
      finished=0;  /* redundant */
      while ((i<sym_last) && (finished==0)) {
	if (symtable[i].tag==spec) {
	  if ((strcmp(name, symtable[i].d.s.name)==0) && 
	      (symtable[i].d.s.charge==charge)) 
	    finished=1;
	}; /* if */
	i++;
      } /* while */
      if (finished==0) {
	symtable[sym_last].tag=spec;
	(void) strcpy(symtable[sym_last].d.s.name, name);
	symtable[sym_last].d.s.charge=charge;
	symtable[sym_last].d.s.conc0=Conc0Default;
	for(j=0; j<MaxSpecConst; j++)
	  symtable[i].d.s.speck[j].used=0;
	symtable[sym_last].d.s.is_param=1;
	sym_last++;
      } /* if */
      finished=1;
      error=NoError;
    } /* if */
  } /* else */
} /* NewParamConc */

int NumOfParameter(void) {

  error=NoError;  
  return (param_last);
} /* NumOfParameter */

void GetParamNo(int no, char *name, double *charge, int *form) {

  if (no>param_last) 
    error=NotFound;
  else {
    strcpy(name, paramtable[no-1].name);
    *charge=paramtable[no-1].charge;
    *form=paramtable[no-1].in_use;
    error=NoError;
  } /* else */
} /* GetParamNo */  

int NumOfConstants(void) {

  int i, no;
  
  no=0;
  for(i=0; i<sym_last; i++)
    if (symtable[i].tag==konst)
      no++;
  error=NoError;
  return no;
} /* NumOfConstants */

void GetConstantNo(int no, char *name) {

  int i;

  i=0;
  while ((i<sym_last) && (no>0)) { 
    if (symtable[i].tag==konst) 
      no--;
    i++;
  } /* while */
  if (no==0) {
    strcpy(name, symtable[i-1].d.k.name);
    error=NoError;
  } /* if */
  else
    error=NotFound;
} /* GetConstantNo */

void NewPrintVar(char *name) {

  if (prn_last>=MaxPrint)
    error=TooManyPrn;
  else {
    prntable[prn_last].tag=1;
    strcpy(prntable[prn_last].name, name);
    prn_last++;
    error=NoError;
  } /* else */
} /* NewPrintVar */

void NewPrintConc(char *name, double charge) {

  if (prn_last>=MaxPrint)
    error=TooManyPrn;
  else {
    prntable[prn_last].tag=2;
    strcpy(prntable[prn_last].name, name);
    prntable[prn_last].charge=charge;
    prn_last++;
    error=NoError;
  } /* else */
} /* NewPrintConc */

int NumOfPrint(void) {
 
  error=NoError;
  return prn_last;
} /* NumOfPrint */

void GetPrintNo(int no, char *name, double *charge, int *tag) {

  if (no>=prn_last)
    error=NotFound;
  else {
    *tag=prntable[no-1].tag;
    strcpy(name, prntable[no-1].name);
    *charge=prntable[no-1].charge;            /* may give unusable values */
    error=NoError;
  } /* else */
} /* GetPrintNo */

int IsSpecInPrnList(char *name, double charge, int tag) {

  /* return 0 is species or dyn. var. should printed */

  int i;

  i=0;
  while (i<prn_last) {
    if (prntable[i].tag==tag) {
      if (tag==1) {             /* dyn var */
	if (strcmp(prntable[i].name, name)==0)
	  return 0;
      } else {         /* tag == 2 i.e. species */
	if ((strcmp(prntable[i].name, name)==0) && 
	    (prntable[i].charge==charge))
	  return 0;
      } /* else */
    } /* if */
    i++;
  } /* while */
  return 1;
} /* IsSpecInPrnList */
	    
int IsSpecParam(char *name, double charge) {

  int i, finished;

  finished=0;
  i=0;
  while ((i<param_last) && (finished==0)) {
    if (paramtable[i].in_use==2) 
      if ((strcmp(paramtable[i].name, name)==0) && 
	  (paramtable[i].charge==charge))
	finished=1;
    i++;
  } /* while */
  return finished;
} /* IsSpecParam */

void NewLowHighPrefParam(char *name, double low, double high, double pref) {

  int i, finished;

  finished=0;
  i=0;
  while ((i<param_last) && (finished==0)) {
    if (paramtable[i].in_use==1) {
      if (strcmp(name, paramtable[i].name)==0) {
	paramtable[i].low=low;
	paramtable[i].high=high;
	paramtable[i].pref=pref;
	finished=1;
	error=NoError;
      } /* if */
    } /* if */
    i++;
  } /* while */
  if (finished==0)
    error=NotFound;
} /* NewLowHighPrefParam */

void NewLowHighPrefConc(char *name, double charge, double low, 
			double high, double pref) {

  int i, finished;

  finished=0;
  i=0;
  while ((i<param_last) && (finished==0)) {
    if (paramtable[i].in_use==2) {
      if ((strcmp(name, paramtable[i].name)==0) && 
	  (charge==paramtable[i].charge)) {
        paramtable[i].low=low;
        paramtable[i].high=high;
        paramtable[i].pref=pref;
        finished=1;
        error=NoError;
      } /* if */
    } /* if */
    i++;
  } /* while */
  if (finished==0)
    error=NotFound;
} /* NewLowHighPrefConc */

void GetLowHighPrefParam(char *name, double *low, double *high, 
			 double *pref) {

  int i, finished;

  finished=0;
  i=0;
  while ((i<param_last) && (finished==0)) {
    if (paramtable[i].in_use==1) {
      if (strcmp(name, paramtable[i].name)==0) {
	*low=paramtable[i].low;
	*high=paramtable[i].high;
	*pref=paramtable[i].pref;
	finished=1;
	error=NoError;
      } /* if */
    }
    i++;
  } /* while */
  if (finished==0)
    error=NotFound;
} /* GetLowHighPrefParam */

void GetLowHighPrefConc(char *name, double charge, double *low, 
			double *high, double *pref) {

  int i, finished;

  finished=0;
  i=0;
  while ((i<param_last) && (finished==0)) {
    if (paramtable[i].in_use==2) {
      if ((strcmp(name, paramtable[i].name)==0) && 
	  (charge==paramtable[i].charge)) {
        *low=paramtable[i].low;
        *high=paramtable[i].high;
        *pref=paramtable[i].pref;
        finished=1;
        error=NoError;
      } /* if */
    } /* if */
    i++;
  } /* while */
  if (finished==0)
    error=NotFound;
} /* GetLowHighConc */
    
void GetInitParam(char *name, double *val) {

  int i, finished;

  finished=0;
  i=0;
  while ((i<param_last) && (finished==0)) {
    if (paramtable[i].in_use==1) {
      if (strcmp(name, paramtable[i].name)==0) {
	*val=paramtable[i].init_value;
	finished=1;
	error=NoError;
      } /* if */
    } /* if */
    i++;
  } /* while */
  if (finished==0)
    error=NotFound;
} /* GetInitParam */

void GetInitConc(char *name, double charge, double *val) {

  int i, finished;

  finished=0;
  i=0;
  while ((i<param_last) && (finished==0)) {
    if (paramtable[i].in_use==2) {
      if ((strcmp(name, paramtable[i].name)==0) && 
	  (charge==paramtable[i].charge)) {
	*val=paramtable[i].init_value;
	finished=1;
	error=NoError;
      } /* if */
    } /* if */
    i++;
  } /* while */
  if (finished==0)
    error=NotFound;
} /* GetInitConc */

void GetDeltaParam(char *name, double *val) {

  int i, finished;

  finished=0;
  i=0;
  while ((i<param_last) && (finished==0)) {
    if (paramtable[i].in_use==1) {
      if (strcmp(name, paramtable[i].name)==0) {
	*val=paramtable[i].delta;
	finished=1;
	error=NoError;
      } /* if */
    } /* if */
    i++;
  } /* while */
  if (finished==0)
    error=NotFound;
} /* GetDeltaParam */

void GetDeltaConc(char *name, double charge, double *val) {

  int i, finished;

  finished=0;
  i=0;
  while ((i<param_last) && (finished==0)) {
    if (paramtable[i].in_use==2) {
      if ((strcmp(name, paramtable[i].name)==0) && 
	  (charge==paramtable[i].charge)) {
        *val=paramtable[i].delta;
        finished=1;
        error=NoError;
      } /* if */
    } /* if */
    i++;
  } /* while */
  if (finished==0)
    error=NotFound;
} /* GetDeltaParam */

void NewDirectForParam(char *name, int val) {

  int i, finished;
  
  finished=0;
  i=0;
  error=NotFound;
  while ((i<param_last) && (finished==0)) {
    if (paramtable[i].in_use==1) {
      if (strcmp(name, paramtable[i].name)==0) {
	paramtable[i].direct=val;
	finished=1;
	error=NoError;
      } /* if name==paramtable.name */
    } /* if in use */
    i++;
  } /* while */
} /* NewDirectForParam */

int GetDirectForParam(char *name) {

  int i, finished;

  finished=0;
  i=0;
  while ((i<param_last) && (finished==0)) {
    if (paramtable[i].in_use==1) {
      if (strcmp(name, paramtable[i].name)==0) {
	error=NoError;
	finished=1;
	return paramtable[i].direct;
      } /* if name == param.name */
    } /* if in use */
    i++;
  } /* while */
  error=NotFound;
  return 0;
} /* GetDirectForParam */

void NewDirectForConc(char *name, double charge, int val) {

  int i, finished;

  finished=0;
  i=0;
  error=NotFound;
  while ((i<param_last) && (finished==0)) {
    if (paramtable[i].in_use==2) {
      if ((strcmp(name, paramtable[i].name)==0) && 
	  (charge==paramtable[i].charge)) {
	finished=1;
	paramtable[i].direct=val;
	error=NoError;
      } /* if name == param.name */
    } /* if in use */
    i++;
  } /* while */
} /* NewDirectForConc */

int GetDirectForConc(char *name, double charge) {
  
  int i, finished;

  finished=0;
  i=0;
  error=NotFound;
  while ((i<param_last) && (finished==0)) {
    if (paramtable[i].in_use==2) {
      if ((strcmp(name, paramtable[i].name)==0) && 
	  (charge==paramtable[i].charge)) {
	finished=1;
	error=NoError;
	return paramtable[i].direct;
      }
    }
    i++;
  }
  return 0;
} /* GetDirectForConc */

void NewInitValue(char *name, double val) {

  int i, found;

  i=0;
  found=0;
  while ((i<expr_last) && (found==0)) {
    if (strcmp(exprtable[i].name, name)==0) {
      error=ExprAlready;
      found=1;
    } /* if */
    else
      i++;
  } /* while */
  if (found==1) {
    exprtable[i].init_value=val;
    error=NoError;
  } /* if */
  else 
    error=NotFound;
} /* NewInitValue */

double GetInitValue(char *name) {

  int i, found;
  double temp;

  error=NotFound;
  i=0;
  found=0;
  while ((i<expr_last) && (found==0)) {
    if (strcmp(exprtable[i].name, name)==0) {
      temp=exprtable[i].init_value;
      found=1;
      error=NoError;
    } /* if */
    else
      i++;
  } /* while */
  return temp;
} /* GetInitValue */

int NoSpecInReacs(char *name, double charge) {

  int i, j, j_finished;
  int res=0;

  i=0;
  while (i<react_last) {
    j_finished=0;
    j=0;
    while ((j<MaxSpec) && (j_finished==0)) {
      if (reacttable[i].species[j].in_use==1)
	if ((strcmp(reacttable[i].species[j].name, name)==0) && 
	    (reacttable[i].species[j].charge==charge)) {
	  res++;
	  j_finished=0;
	}
      j++;
    }
    i++;
  }
  return res;
} /* NoSpecInReacs */

void NewStrConst(char *name, char *value) {

  if (str_last>=StrConstSize)
    error=TooManyConst;
  else {
    strcpy(strtable[str_last].name, name);
    strcpy(strtable[str_last].value, value);
    str_last++;
    error=NoError;
  }
} /* NewStrConst */

void GetStrConst(char *name, char *value) {

  int i, finished;

  i=0;
  finished=0;
  while ((i<str_last) && (finished==0)) {
    if (strcmp(name, strtable[i].name)==0) {
      finished=1;
      error=NoError;
      strcpy(value, strtable[i].value);
    }
    i++;
  }
  if (finished==0)
    error=NotFound;
} /* GetStrConst */ 

int NumOfStrConst(void) {

  return str_last;
} /* NumOfStrConst */

void GetStrConstNo(int no, char *name, char *value) {

  if (no>str_last)
    error=NotFound;
  else {
    strcpy(value, strtable[no-1].value);
    strcpy(name, strtable[no-1].name);
    error=NoError;
  }
} /* GetStrConstNo */

void NonAutoSystem(void) {

  autonom_system=1;
} /* NonAutoSystem */

int IsNonAutoSystem(void) {

  return autonom_system;
} /* IsNonAutoSystem */

void GetStocMatrix(int n, int m, double **nu) {

  int     i, j, reactno;
  double  charge, coeff;
  char    name[STRING_LENGTH];
  

  for(i=0; i<n; i++) {
    reactno=GetReactNo(i);
    for(j=0; j<m; j++) {
      GetSpecNo(i, name, &charge);
      coeff=GetCoeffInReact(reactno, name, charge, 1);
      coeff-=GetCoeffInReact(reactno, name, charge, 2);
      nu[i][j]=coeff;
    } /* for j */
  } /* for i */
} /* GetStocMatrix */
      
