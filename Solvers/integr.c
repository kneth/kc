/****************************************************************************
  This is a small library for integrating functions.

  All the integration routines assume that the integrand has the form
     double f(double x, double *params)
  where x is the independent variable and params is a vector of parameters.

  CopyWrong 1994 by
  Kenneth Geisshirt (kneth@osc.kiku.dk)
  Department of Theoretical Chemistry
  H.C. Orsted Institute
  Universitetsparken 5
  2100 Copenhagen
  Denmark

  References:
  [1]  Numerical Analysis.
       D. Kincaid, and W. Cheney, Brooks/Cole, 1991.


  Last updated: 30 August 1994
******************************************************************************/

#include "integr.h"
#include <math.h>

/*****************************************************************************
  Gauss5 is a five-point Gauss quadrature, [1, pp.456-462]. The argument n
  is the number of subintervals.
******************************************************************************/

double Gauss5(int n, double a, double b, double *params, 
	      double (*f)(double, double *)) {

  const double x_0 = 0.0;
  const double x_1 = 0.538469310105683;
  const double x_2 = 0.906179845938664;
  const double w_0 = 0.568888888888889;
  const double w_1 = 0.478628670499366;
  const double w_2 = 0.236926885056189;

  double u, v, S, inte, h;
  int    i;

  h=(b-a)/((double) n);
  inte=0.0;

  for(i=0; i<=(n-1); i++) {
    u=(h*x_0+((double) 2*i+1)*h)*0.5;
    S=w_0*f(u, params);
    u=(h*x_1+((double) 2*i+1)*h)*0.5;
    v=(-h*x_1+((double) 2*i+1)*h)*0.5;
    S+=w_1*(f(u, params)+f(v, params));
    u=(h*x_2+((double) 2*i+1)*h)*0.5;
    v=(-h*x_2+((double) 2*i+1)*h)*0.5;
    S+=w_2*(f(u, params)+f(v, params));
    S*=h*0.5;
    inte+=S;
  } /* for i */
  return inte;
} /* Gauss5 */


/*****************************************************************************
  Trapez is using the trapeziod rule in order to evaluate the integral, 
  [1, pp. 445-446].

  The parameter n is the number of subintervals.
******************************************************************************/
  
double Trapez(int n, double a, double b, double *params, 
	      double (*f)(double, double *)) {

  int    i;         /* counter               */
  double h;         /* length of subinterval */
  double sum;       /* temporary real        */
  double x;         /* counter               */

  h=(b-a)/((double) n);          /* find subintervals */
  sum=f(a, params)+f(b, params); /* init.             */
  x=a+h;
  for(i=1; i<=(n-1); i++) {
    sum+=2.0*f(x, params);
    x+=h;
  } /* for i */
  sum*=0.5*h;
  return sum;
} /* Trapez */


/*****************************************************************************
  Simpson is computing the integral by applying composite Simpson's rule, 
  [1, pp. 447-448].
******************************************************************************/

double Simpson(int n, double a, double b, double *params, 
	       double (*f)(double, double *)) {

  int         i;   /* counter                */
  double      h;   /* length of subintervals */
  double      sum; /* sum so far             */
  
  h=(b-a)/((double) n);
  sum=f(a, params)+f(b, params);
  for(i=2; i<=n/2; i++) 
    sum+=2.0*f(a+((double)2*i-2)*h, params);
  for(i=1; i<=n/2; i++)
    sum+=4.0*f(a+((double)2*i-1)*h, params);
  sum*=h/3.0;
  return sum;
} /* Simpson */
