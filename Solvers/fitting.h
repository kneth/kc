/**************************************************************************
  Header file for the fitting routines in SSL.

  (C) Copyright 1996 by
  Kenneth Geisshirt (kneth@fatou.ruc.dk)      Keld Nielsen (kn@kin.kiku.dk)
  Dept. of Life Sciences and Chemistry       Dept. of Theoretical Chemistry
  Roskilde University                              University of Copenhagen
  Marbjergvej 35                                       Universitetsparken 5
  DK-4000 Roskilde                                       DK-2100 Copenhagen

  Last modified: 24 March 1996 by KG
***************************************************************************/

#ifndef _FIT_LIB_
#define _FIT_LIB_

extern void LineFit(const int, double *, double *, double *, double *, 
                    double *, double *);'
extern void LMfit(const int, const int,   );

#endif /* _FIT_LIB_ */
