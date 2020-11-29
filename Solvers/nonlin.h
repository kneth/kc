/***************************************************************************
  (C) Copyright 1995 by
  Kenneth Geisshirt (kneth@fatou.ruc.dk)
  Department of Life Sciences and Chemistry
  Roskilde University
  P.O. Box 260
  4000 Roskilde
  Denmark

  Voice: +45 46 75 77 11, ext. 27 45
  Fax:   +45 46 75 77 21
  

  This is the header file for a library containing solvers of nonlinear
  algebraic equations.

  Last modified: 27th September 1995
*****************************************************************************/

#ifndef _NONLIN_LIB_
#define _NONLIN_LIB_

extern double NewtonRaphson1D(double, double (*f)(double), 
			      double (*df)(double), double, int);
extern void NewtonRaphson(int, double *, double *, 
			  void (*f)(double *, double *), 
			  void (*df)(double *, double **), double, int);

#endif
