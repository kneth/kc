/****************************************************************************
  CodeGen2 - new code generator routines for kc.
  
  (C) Copyright 1992-1996 by 
  Kenneth Geisshirt (kneth@fatou.ruc.dk)        Keld Nielsen (kn@kin.kiku.dk)
  Dept. of Life Sciences and Chemistry         Dept. of Theoretical Chemistry
  Roskilde University                                University of Copenhagen
  P.O. Box 260                                           Universitetsparken 5
  DK-4000 Roskilde                                    DK-2100 Copenhagen East

  See kc.tex for details.

  Last updated: 24 May 1996 by KG
*****************************************************************************/


#include "codegen2.h"
#include "symbmath.h"
#include "tableman.h"
#include "misc.h"

#include <stdlib.h>
#include <stdio.h>

/****************************************************************************
  GenRateAndJac calculates the rate expressions and the jacobian matrix. The
  storages needed is also allocated by this routine and the communication
  is done by global variables.
*****************************************************************************/

void GenRateAndJac() {

  
  int    numcon;             /* number of constaints        */
  int    numspec;            /* number of species           */
  int    numdyn;             /* number of dyn. variables    */
  int    numreact;           /* number of reactions         */

  Tree   tmp;                /* temporary expression        */
  char   *name;              /* species/dyn.var. name       */
  double charge;             /* charge of species           */



  /***** A few definitions                              *****/
  numcon=NumOfConstraint();
  numspec=NoOfSpec();
  numdyn=NoOfDynVar();
  numreact=NoOfReact();

  
  /***** Allocate storage                               *****/
  name=StringAlloc();

  stoccon=(Tree *)calloc(numcon, sizeof(Tree));
  if (stoccon==NULL) fprintf(stderr, "GenRateAndJac: Not enough space\n");

  rate=(Tree *)calloc(numdyn+numreact, sizeof(Tree));
  if (rate==NULL) fprintf(stderr, "GenRateAndJac: Not enough space\n");

  stocdiff=(Tree **)calloc(numcon, sizeof(Tree *));
  if (stocdiff==NULL) fprintf(stderr, "GenRateAndJac: Not enough space\n");
  for(i=0; i<numcon; i++) {
    stocdiff[i]=(Tree *)calloc(numspec+numdyn, sizeof(Tree));
    if (stocdiff[i]==NULL) 
      fprintf(stderr, "GenRateAndJac: Not enough space\n");
  } /* for i=0..numcon-1 */

  stomatrix=(double **)calloc(numspec+numdyn, sizeof(Tree *));
  if (stomatrix==NULL) fprintf(stderr, "GenRateAndJac: Not enough space\n");
  for(i=0; i<(numspec+numdyn); i++) {
    stomatrix[i]=(double *)calloc(numreact, sizeof(double));
    if (stomatrix[i]==NULL) 
      fprintf(stderr, "GenRateAndJac: Not enough space\n");
  } /* for i=0..numspec+numdyn-1 */


  /***** Process constaints                             *****/
  for(i=0; i<numcon; i++) { 
    tmp=TreeCreate();
    stoccon[i]=TreeCreate();
    GetConstraintNo(i+1, name, &charge, tmp);
    RenameSpec(rename, name, charge);
    for(j=0; j<i; j++) 
      TreeSubstTree(con[j], rename, tmp);
    TreeCpy(stoccon[i], tmp);
    TreeKill(tmp);
  } /* for i=0..numcon-1 */


  /***** Compute the rate expressions                   *****/
  
