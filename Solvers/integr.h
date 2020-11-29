/****************************************************************************
  This is a small library for integrating functions.
  Documentation and references are found in integr.c.

  CopyWrong 1994 by
  Kenneth Geisshirt (kneth@osc.kiku.dk)
  Department of Theoretical Chemistry
  H.C. Orsted Institute
  Universitetsparken 5
  2100 Copenhagen
  Denmark

  Last updated: 30 August 1994
******************************************************************************/

#ifndef _INTEGR_LIB_
#define _INTEGR_LIB_

extern double   Gauss5(int, double, double, double *, 
		       double (*f)(double, double *));
extern double   Trapez(int, double, double, double *, 
		       double (*f)(double, double *));
extern double   Simpson(int, double, double, double *, 
			double (*f)(double, double *));
#endif
