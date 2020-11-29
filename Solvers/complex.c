/*****************************************************************************
  Complex - a small library for doing complex algebra in C.

  CopyWrong 1994 by
  Kenneth Geisshirt (kneth@osc.kiku.dk)
  Department of Theoretical Chemstry
  H.C. Orsted Institute
  Universitetsparken 5
  2100 Copenhagen 5
  Denmark

  Last updated: 12 September 1994
******************************************************************************/


#include <math.h>
#include "complex.h"

void ComplexAssign(double re, double im, Complex *z) {

  z->re=re;
  z->im=im;
} /* ComplexAssign */ 
  

void ComplexAdd(Complex z1, Complex z2, Complex *res) {

  res->re=z1.re+z2.re;
  res->im=z1.im+z2.im;
} /* ComplexAdd */


void ComplexSub(Complex z1, Complex z2, Complex *res) {

  res->re=z1.re-z2.re;
  res->im=z1.im-z2.im;
} /* ComplexSub */
 

void ComplexMul(Complex z1, Complex z2, Complex *res) {

  res->re=z1.re*z2.re-z1.im*z2.im;
  res->im=z1.re*z2.im+z1.im*z2.re;
} /* ComplexMul */


void ComplexDiv(Complex a, Complex b, Complex *res) {

  double temp;

  temp=b.re*b.re+b.im*b.im;
  res->re=(a.re*b.re+a.im*b.im)/temp;
  res->im=(a.im*b.re-a.re*b.im)/temp;
} /* ComplexDiv */


double ComplexNorm(Complex z) {

  return (sqrt(z.re*z.re+z.im*z.im));
} /* ComplexNorm */


double ComplexArg(Complex z) {

  return (atan2(z.im, z.re));
} /* ComplexArg */ 

