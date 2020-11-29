/*****************************************************************************
  Header file for service routines for ODE solvers.

  CopyWrong 1994-1995 by
  Kenneth Geisshirt (kneth@fatou.ruc.dk) 
  Department of Life Sciences and Chemistry
  Roskilde University
  P.O. Box 260
  4000 Roskilde
  Denmark
  
  and
  Keld Nielsen (kn@osc.kiku.dk)
  Department of Theoretical Chemistry
  H.C. Orsted Institute
  Universitetsparken 5
  2100 Copenhagen
  Denmark

  Last updated: 15 February 1995 by KG
******************************************************************************/

#ifndef _ODE_SERV_LIB_
#define _ODE_SERV_LIB_

#include <stdio.h>

extern void     kcerror(char *);
extern double   FindMachEps(void);
extern double   MaxVec(int, double *);
extern double   MaxPair(double, double);
extern double   MinPair(double, double);
extern void     PrintState(int, int, int, double, double *, FILE *, int *);
extern void     BSort(const int, double *);
extern double   KronDelta(int, int);
#endif


