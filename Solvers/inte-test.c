/***************************************************************************
  Just a small test of the integration library.

  Kenneth Geisshirt, 30 August 1994.
****************************************************************************/

#include <math.h>
#include <stdio.h>
#include "integr.h"

double f(double x, double params[0]) {

  return sin(x);
}

void main(void) {

  double param[0];
  
  printf("Trapez  (50): %2.10e\n", Trapez(50, 0.0, PI*0.5, param, &f));
  printf("Simpson (50): %2.10e\n", Simpson(50, 0.0, PI*0.5, param, &f));
  printf("Gauss5  (10): %2.10e\n", Gauss5(10, 0.0, PI*0.5, param, &f));
} 
