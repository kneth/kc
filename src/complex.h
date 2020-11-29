/*****************************************************************************
  Complex - a small library for doing complex algebra in C.

  CopyWrong 1994 by
  Kenneth Geisshirt (kneth@osc.kiku.dk)
  Department of Theoretical Chemstry
  H.C. Orsted Institute
  Universitetsparken 5
  2100 Copenhagen 5
  Denmark

  Last updated: 20 August 1994
******************************************************************************/

#ifndef _COMPLEX_LIB_
#define _COMPLEX_LIB_

struct ComplexStruct {
  double   re, im;
};

typedef struct ComplexStruct Complex;

extern void    ComplexAssign(double, double, Complex *);
extern void    ComplexAdd(Complex, Complex, Complex *);
extern void    ComplexSub(Complex, Complex, Complex *); 
extern void    ComplexMul(Complex, Complex, Complex *);
extern void    ComplexDiv(Complex, Complex, Complex *);
extern double  ComplexNorm(Complex);
extern double  ComplexArg(Complex);

#endif
