/************************************************************************
  A code generator which is using VODE [1] as integrator.

  CopyWrong 1994 by
  Kenneth Geisshirt (kneth@osc.kiku.dk)
  Department of Theoretical Chemistry
  H.C. Orsted Institute
  Universitetsparken 5
  2100 Copenhagen
  Denmark

  Reference:
  [1] VODE: A Variable-Coefficient ODE Solver.
      P.N. Brown et al., SIAM J. Sci. Stat. Comp.
      pp. 1038-1051, volume 10, number 5, 1989.

  Last updated: 23 April 1994
************************************************************************/

#include <stdio.h>

extern void VODE(FILE *);
