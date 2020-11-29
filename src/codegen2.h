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


#ifndef _CODE_GEN2_H_
#define _CODE_GEN2_H_

/***** Global variables                               *****/
Tree   *stoccon;           /* stoichiometric constraints  */
Tree   *rate;              /* rate expressions            */
Tree   **stocdiff;         /* diff. of constaints         */
double **stomatrix;        /* stoichiometric matrix       */



#endif
