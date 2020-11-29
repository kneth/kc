/*****************************************************************************
  Quench is a library for calculation of quenching data.
  See quench.h for details.
  Last updated: 16 May 1995  by KN
*****************************************************************************/

#ifndef _QUENCH_LIB_
#define _QUENCH_LIB_

#include "complex.h"
#include "matrix.h"
#include <math.h>
#include <stdio.h>

extern double arctan_local(Complex);
extern void compamppha(int, double **, double *, double *);
extern void stopdata(int, int, double **, double *,
	      double *, double *, double *, double *);

#endif
