/*****************************************************************************
  Header file for ODE solver library.

  CopyWrong 1994 by
  Kenneth Geisshirt (kneth@osc.kiku.dk) and
  Keld Nielsen (kn@osc.kiku.dk)
  Department of Theoretical Chemistry
  H.C. Orsted Institute
  Universitetsparken 5
  2100 Copenhagen
  Denmark

  Last updated: 6 September 1994
******************************************************************************/

#ifndef _ODE_SOLV_LIB_
#define _ODE_SOLV_LIB_

double **jacobi_matx;               /* global variable */

extern void  CalahanOneStep(int, double, double *, double *, 
			    void (*f_)(double *, double *), 
			    void (*jac_)(double *));
extern void  Calahan(int, double, double *, double *, double *, 
		     void (*f_)(double *, double *), void (*jac_)(double *));
extern void  RKFNC(int, double, double *, double *, double *, 
		   void (*f_)(double *, double *));
extern void  RK4OneStep(int, double, double *, double *, 
			void (*f_)(double *, double *));
extern void  RK4(int, double, double *, double *, double *,
		 void (*f_)(double *, double *));
extern void  GRK4T(int, double, double *, double *, double *, 
		   void (*f_)(double *, double *), void (*jac_)(double *));
extern void  RKFNCTime(int, int, double, double, double *, double *, double *,
		       void (*f_)(double *, double *));
extern void  GRK4TTime(int, int, double, double, double *, double *, double *,
		       void (*f_)(double *, double *), void (*jac_)(double *));

#endif
