/***************************************************************************
  This source file contains various numerical solvers for ordinary
  differential equations.

  The routines associated by each solver compute one and only one time
  step.

  Please see the file kksolver.c for all the details. The file is also
  an example of how this library can be used.


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


  References:
  [1]  W.H. Press, et al. Numerical Recipes in C, 2. edition.
       Cambridge  University Press, 1992.
  [2]  D.A. Calahan. (1968). Proc. IEEE, April 1968, p. 744.
  [3]  P. Kaps, and P. Rentrop. (1979). Numer. Math, 33, pp. 55-68.
  [4]  K. Geisshirt. Chemical Waves in Reaction-Diffusion Systems: A Numerical 
       Study. (M.Sc. thesis, University of Copenhagen), 1994.
  [5]  M. Kubicek, and M. Marek, Computational Methods In Bifucation Theory
       And Dissipative Structures, Springer-Verlag, New York, 1983.
  [6]  J.R. Cash, and A.H. Karp. (1990). ACM Trans. Math. Softw. 16, 
       pp. 201-222.
  [7]  P. Kaps, S.W.H. Poon, and T.D. Bui. (1985). Computing. 34, pp. 17-40.
  [8]  P. Deuflhard, G. Bader, and U. Nowak. (1981). pp. 38-55.
       In K.H. Ebert, P. Deuflhard, and W. Jager (eds.). Modelling of Chemical
       Reaction Systems. Springer Series in Chemical Physics, Vol. 18.
       Springer-Verlag, Berlin, 1981.
  [9]  E. Hairer, and G. Wanner. Solving Ordinary Differential Equations II.
       Springer Series In Computational Mathematics, Vol. 14. Springer-Verlag,
       Berlin, 1991.
  [10] D. Kincaid, and W. Cheney. Numerical Analysis. Brooks/Cole, 1991.

  Last updated: 15 February 1995 by KG
****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <time.h>
#include "odesolv.h"

/***************************************************************************
  The matrix manipulating routines like linear equation solvers are all
  implemented in a small library. It is imported below. 
****************************************************************************/

#include "matrix.h"

/****************************************************************************
  Solver_ = 1:  For stiff autonomous systems, [2, 5].

  Semi-implicit third order Rosenbrock method. The constants a1_, b1_, R1_,
  and R2_ are due to Calahan (e.g. [2] and [5]).
  Adaptive stepsize control: Step-doubling.
  The routine CalahanOneStep take one step of the integration algorithm, and
  the routine Calahan take care of the step-doubling.
  The output from Calahan: xnew_ and xerr_ contain respectively the new values
  and the estimate of the local error of the dependent variables.
*****************************************************************************/

void CalahanOneStep(int N_, double htime_, double *x_, double *xnew_, 
		    void (*f_)(double *, double *), void (*jac_)(double *)) {

  static const double a1_ = 0.788675134594813; /* (3.0+sqrt(3.0))/6.0 */
  static const double b1_ = -1.15470053837925; /* -2.0/sqrt(3.0) */
  static const double R1_ = 0.75;
  static const double R2_ = 0.25;

  double **A_;
  double *k1_, *k2_, *temp_;
  double dummy_;
  int    *index;
  int    i, j;

  A_=MatrixAlloc(N_);
  index=IntVectorAlloc(N_);
  temp_=VectorAlloc(N_);
  k1_=VectorAlloc(N_);
  k2_=VectorAlloc(N_);

  dummy_=1.0/htime_;
  f_(x_, k1_);
  jac_(x_);
  for(i=0; i<N_; i++)
    for(j=0; j<N_; j++)
      A_[i][j]=-a1_*jacobi_matx[i][j];
  for(i=0; i<N_; i++)
    A_[i][i]=A_[i][i]+dummy_;
  LUfact(N_, A_, index);
  LUsubst(N_, A_, index, k1_);
  for(i=0; i<N_; i++)
    temp_[i]=x_[i]+b1_*k1_[i];
  f_(temp_, k2_);
  LUsubst(N_, A_, index, k2_);
  for(i=0; i<N_; i++)
    xnew_[i]=x_[i]+R1_*k1_[i]+R2_*k2_[i];

  MatrixFree(N_, A_);
  IntVectorFree(N_, index);
  VectorFree(N_, temp_);
  VectorFree(N_, k1_);
  VectorFree(N_, k2_);
} /* CalahanOneStep */

void Calahan(int N_, double htime_, double *x_, double *xnew_, double *xerr_, 
	     void (*f_)(double *, double *), void (*jac_)(double *)) {

  double *xh_, *xh2_, halfhtime_;
  int    i;

  xh_=VectorAlloc(N_);
  xh2_=VectorAlloc(N_);

  halfhtime_=0.5*htime_;
  CalahanOneStep(N_, halfhtime_, x_, xh_, f_, jac_);
  CalahanOneStep(N_, halfhtime_, xh_, xh2_, f_, jac_);
  CalahanOneStep(N_, htime_, x_, xnew_, f_, jac_);
  for(i=0; i<N_; i++)
    xerr_[i]=(xnew_[i]-xh2_[i]);

  VectorFree(N_, xh_);
  VectorFree(N_, xh2_);
} /* Calahan */


/****************************************************************************
  Solver_ = 2: (see also solver_ 5):  For non-stiff autonomous systems.
  ([6] and [1, pp. 711-724])

  A fixed-fifth order Runge-Kutta method. The constants are due to Cash
  and Karp ([6] and [1, pp. 711-724]).
  Adaptive stepsize control: Estimate of local error due to embedded fourth
  order method.
  The output from RKFNK: xnew_ and xerr_ contain respectively the new values
  and the estimate of the local error of the dependent variables.
*****************************************************************************/

void RKFNC(int N_, double htime_, double *x_, double *xnew_, double *xerr_, 
	   void (*f_)(double *, double *)) {

  double *k1_, *k2_, *k3_, *k4_, *k5_, *k6_;

  static const double   b21_ =      1.0/5.0;
  static const double   b31_ =      3.0/40.0;
  static const double   b32_ =      9.0/40.0;
  static const double   b41_ =      3.0/10.0;
  static const double   b42_ =     -9.0/10.0;
  static const double   b43_ =      6.0/5.0;
  static const double   b51_ =    -11.0/54.0;
  static const double   b52_ =      5.0/2.0;
  static const double   b53_ =    -70.0/27.0;
  static const double   b54_ =     35.0/27.0;
  static const double   b61_ =   1631.0/55296.0;
  static const double   b62_ =    175.0/512.0;
  static const double   b63_ =    575.0/13824.0;
  static const double   b64_ =  44275.0/110592.0;
  static const double   b65_ =    253.0/4096.0;
  static const double   c1_  =     37.0/378.0;
  static const double   dc1_ =     37.0/378.0 -  2825.0/27648.0;
  static const double   c3_  =    250.0/621.0;
  static const double   dc3_ =    250.0/621.0 - 18575.0/48384.0;
  static const double   c4_  =    125.0/594.0;
  static const double   dc4_ =    125.0/594.0 - 13525.0/55296.0;
  static const double   dc5_ =      0.0 - 277.0/14336.0;
  static const double   c6_  =    512.0/1771.0;
  static const double   dc6_ =    512.0/1771.0 - 1.0/4.0;

  double *temp1_;
  int    i;

  k1_=VectorAlloc(N_);
  k2_=VectorAlloc(N_);
  k3_=VectorAlloc(N_);
  k4_=VectorAlloc(N_);
  k5_=VectorAlloc(N_);
  k6_=VectorAlloc(N_);
  temp1_=VectorAlloc(N_);

  f_(x_, k1_);
  for(i=0; i<N_; i++) {
    k1_[i]*=htime_;
    temp1_[i]=x_[i]+b21_*k1_[i];
  } /* for i */
  f_(temp1_, k2_);
  for(i=0; i<N_; i++) {
    k2_[i]*=htime_;
    temp1_[i]=x_[i]+b31_*k1_[i]+b32_*k2_[i];
  } /* for i */
  f_(temp1_, k3_);
  for(i=0; i<N_; i++) {
    k3_[i]*=htime_;
    temp1_[i]=x_[i]+b41_*k1_[i]+b42_*k2_[i]+b43_*k3_[i];
  } /* for i */
  f_(temp1_, k4_);
  for(i=0; i<N_; i++) {
    k4_[i]*=htime_;
    temp1_[i]=x_[i]+b51_*k1_[i]+b52_*k2_[i]+b53_*k3_[i]+b54_*k4_[i];
  } /* for i */
  f_(temp1_, k5_);
  for(i=0; i<N_; i++) {
    k5_[i]*=htime_;
    temp1_[i]=x_[i]+b61_*k1_[i]+b62_*k2_[i]+b63_*k3_[i]
      +b64_*k4_[i]+b65_*k5_[i];
  } /* for i */
  f_(temp1_, k6_);
  for(i=0; i<N_; i++) {
    k6_[i]*=htime_;
    xnew_[i]=x_[i]+c1_*k1_[i]+c3_*k3_[i]+c4_*k4_[i]+c6_*k6_[i];
    xerr_[i]=dc1_*k1_[i]+dc3_*k3_[i]+dc4_*k4_[i]+dc5_*k5_[i]+dc6_*k6_[i];
  } /* for i */

  VectorFree(N_, k1_);
  VectorFree(N_, k2_);
  VectorFree(N_, k3_);
  VectorFree(N_, k4_);
  VectorFree(N_, k5_);
  VectorFree(N_, k6_);
  VectorFree(N_, temp1_);
} /* RKFNC */


/****************************************************************************
  Solver_ = 3: For non-stiff autonomous systems, e.g. [1, pp. 710-714].

  A fixed-fourth order Runge-Kutta method.
  Adaptive stepsize control: Step-doubling.
  The routine RK4OneStep take one step of the integration algorithm, and
  the routine RK4 take care of the step-doubling.
  The output from RK4: xnew_ and xerr_ contain respectively the new values
  and the estimate of the local error of the dependent variables.
*****************************************************************************/

void RK4OneStep(int N_, double htime_, double *x_, double *xnew_, 
		void (*f_)(double *, double *)) {

  double *k1_, *k2_, *k3_, *k4_;
  double *temp1_, *temp2_;
  int    i;

  k1_=VectorAlloc(N_);
  k2_=VectorAlloc(N_);
  k3_=VectorAlloc(N_);
  k4_=VectorAlloc(N_);
  temp1_=VectorAlloc(N_);
  temp2_=VectorAlloc(N_);

  f_(x_, temp1_);
  for(i=0; i<N_; i++) {
    k1_[i]=htime_*temp1_[i];
    temp1_[i]=x_[i]+0.5*k1_[i];
  } /* for i */
  f_(temp1_, temp2_);
  for(i=0; i<N_; i++) {
    k2_[i]=htime_*temp2_[i];
    temp1_[i]=x_[i]+0.5*k2_[i];
  } /* for i */
  f_(temp1_, temp2_);
  for(i=0; i<N_; i++) {
    k3_[i]=htime_*temp2_[i];
    temp1_[i]=x_[i]+0.5*k3_[i];
  } /* for i */
  f_(temp1_, temp2_);
  for(i=0; i<N_; i++) 
    k4_[i]=htime_*temp2_[i];
  for(i=0; i<N_; i++)
    xnew_[i]=x_[i]+(1.0/6.0)*k1_[i]+(1.0/3.0)*k2_[i]
      +(1.0/3.0)*k3_[i]+(1.0/6.0)*k4_[i];

  VectorFree(N_, k1_);
  VectorFree(N_, k2_);
  VectorFree(N_, k3_);
  VectorFree(N_, k4_);
  VectorFree(N_, temp1_);
  VectorFree(N_, temp2_);
} /* RK4OneStep */

void RK4(int N_, double htime_, double *x_, double *xnew_, 
	 double *xerr_, void (*f_)(double *, double *)) {

  double *xh_, *xh2_, halfhtime_;
  int    i;

  xh_=VectorAlloc(N_);
  xh2_=VectorAlloc(N_);

  halfhtime_=0.5*htime_;
  RK4OneStep(N_, halfhtime_, x_, xh_, f_);
  RK4OneStep(N_, halfhtime_, xh_, xh2_, f_);
  RK4OneStep(N_, htime_, x_, xnew_, f_);
  for(i=0; i<N_; i++)
    xerr_[i]=(xnew_[i]-xh2_[i]);

  VectorFree(N_, xh_);
  VectorFree(N_, xh2_);
} /* RK4 */


/***************************************************************************
  Solver_ = 4. (see also solver_ 7): For stiff autonomous systems, 
  [1, 3, 7, 9].

  Semi-implicit fourth order Rosenbrock method. The constants are due to
  Kaps and Rentrop (version GRK4T ,[3]), but have been transformed by the
  the method outlined in [1, pp. 738-742], [7], and [9].
  Adaptive stepsize control: Estimate of local error due to embedded third
  order method.
  The output from GRK4T: xnew_ and xerr_ contain respectively the new values
  and the estimate of the local error of the dependent variables.
****************************************************************************/

void GRK4T(int N_, double htime_, double *x_, double *xnew_,
	   double *xerr_, void (*f_)(double *, double *),
	   void (*jac_)(double *)) {

  static const double gamma_    =   0.23100000000;
  static const double c21_      =  -5.07167544877;
  static const double c31_      =   6.02015272865;
  static const double c32_      =   0.159750684673;
  static const double c41_      =  -1.856343618677;
  static const double c42_      =  -8.50538085819;
  static const double c43_      =  -2.08407513602;
  static const double a21_      =   2.00000000000;
  static const double a31_      =   4.52470820736;
  static const double a32_      =   4.16352878860;
  /* a41_ = a31_, a42_ = a32_, a43_ = 0.0 */
  static const double m1_       =   3.95750374663;
  static const double m2_       =   4.62489238836;
  static const double m3_       =   0.617477263873;
  static const double m4_       =   1.282612945268;
  static const double dm1_      =  -2.30215540292;
  static const double dm2_      =  -3.07363448539;
  static const double dm3_      =   0.873280801802;
  static const double dm4_      =   1.282612945268;
 
  double **A;
  double *k1_, *k2_, *k3_, *k4_;
  double *temp1_, *temp2_, dummy;
  int    *index;
  int    i, j;

  A=MatrixAlloc(N_);
  index=IntVectorAlloc(N_);
  k1_=VectorAlloc(N_);
  k2_=VectorAlloc(N_);
  k3_=VectorAlloc(N_);
  k4_=VectorAlloc(N_);
  temp1_=VectorAlloc(N_);
  temp2_=VectorAlloc(N_);

  dummy=1.0/(htime_*gamma_);
  jac_(x_);
  for(i=0; i<N_; i++)
    for(j=0; j<N_; j++)
      A[i][j]= -jacobi_matx[i][j];
  for(i=0; i<N_; i++)
    A[i][i]+=dummy;
  LUfact(N_, A, index);

  f_(x_, k1_);
  LUsubst(N_, A, index, k1_);

  for(i=0; i<N_; i++) 
    temp1_[i]=x_[i]+a21_*k1_[i];
  f_(temp1_, k2_);
  for(i=0; i<N_; i++)
    k2_[i]+=((c21_*k1_[i])/htime_);
  LUsubst(N_, A, index, k2_);

  for(i=0; i<N_; i++) 
    temp1_[i]=x_[i]+a31_*k1_[i]+a32_*k2_[i];
  f_(temp1_, k3_);
  for(i=0; i<N_; i++)  {
    k4_[i]=k3_[i];
    k3_[i]+=((c31_*k1_[i]+c32_*k2_[i])/htime_);
  } /* for i */
  LUsubst(N_, A, index, k3_);

  for(i=0; i<N_; i++)
    k4_[i]+=((c41_*k1_[i]+c42_*k2_[i]+c43_*k3_[i])/htime_);
  LUsubst_fixed(k4_);

  for(i=0; i<N_; i++) {
    xnew_[i]= x_[i]+m1_*k1_[i]+m2_*k2_[i]+m3_*k3_[i]+m4_*k4_[i];
    xerr_[i]= dm1_*k1_[i]+dm2_*k2_[i]+dm3_*k3_[i]+dm4_*k4_[i];
  } /* for i */

  MatrixFree(N_, A);
  IntVectorFree(N_, index);
  VectorFree(N_, k1_);
  VectorFree(N_, k2_);
  VectorFree(N_, k3_);
  VectorFree(N_, k4_);
  VectorFree(N_, temp1_);
  VectorFree(N_, temp2_);
} /* GRK4T */

/****************************************************************************
  Solver_ = 5 (see also solver_ 2): For non-stiff non-autonomous systems,
  [6] and [1, pp. 711-724].

  A fixed-fifth order Runge-Kutta method. The constants are due to Cash
  and Karp, [6] and [1, pp. 711-724].
  Adaptive stepsize control: Estimate of local error due to embedded fourth
  order method.
  The output from RKFNKTime: xnew_ and xerr_ contain respectively the new 
  values and the estimate of the local error of the dependent variables.
*****************************************************************************/

void RKFNCTime(int N_, int ns_, double time_, double htime_, double *x_,
	       double *xnew_, double *xerr_, void (*f_)(double *, double *)) {

  double *k1_, *k2_, *k3_, *k4_, *k5_, *k6_;

  /* The constants a2_, a3_, a4_, a5_, and a6_ are only relevant for    */
  /* the Runge-Kutta-Cash-Karp integration methode when the system of   */
  /* differential equations is nonautonomous. In our implementation of  */
  /* the algorithm we do not use these constants, and they are only     */
  /* included for the user that have an interest in changing the source */
  /* code. If one have a nonautonomous system we prefere to transform   */
  /* the system to a system of autonomous differential equations by     */
  /* including a differential equation for the time (the dimension of   */
  /* the system is increased by one): x(n+1) = time, x'(n+1) = 1.0      */

  static const double   a2_  =  1.0/5.0; 
  static const double   a3_  =  3.0/10.0;
  static const double   a4_  =  3.0/5.0;
  static const double   a5_  =  1.0;
  static const double   a6_  =  7.0/8.0;

  static const double   b21_ =      1.0/5.0;
  static const double   b31_ =      3.0/40.0;
  static const double   b32_ =      9.0/40.0;
  static const double   b41_ =      3.0/10.0;
  static const double   b42_ =     -9.0/10.0;
  static const double   b43_ =      6.0/5.0;
  static const double   b51_ =    -11.0/54.0;
  static const double   b52_ =      5.0/2.0;
  static const double   b53_ =    -70.0/27.0;
  static const double   b54_ =     35.0/27.0;
  static const double   b61_ =   1631.0/55296.0;
  static const double   b62_ =    175.0/512.0;
  static const double   b63_ =    575.0/13824.0;
  static const double   b64_ =  44275.0/110592.0;
  static const double   b65_ =    253.0/4096.0;
  static const double   c1_  =     37.0/378.0;
  static const double   dc1_ =     37.0/378.0 -  2825.0/27648.0;
  static const double   c3_  =    250.0/621.0;
  static const double   dc3_ =    250.0/621.0 - 18575.0/48384.0;
  static const double   c4_  =    125.0/594.0;
  static const double   dc4_ =    125.0/594.0 - 13525.0/55296.0;
  static const double   dc5_ =      0.0 - 277.0/14336.0;
  static const double   c6_  =    512.0/1771.0;
  static const double   dc6_ =    512.0/1771.0 - 1.0/4.0;

  double *temp1_;
  int    i;

  k1_=VectorAlloc(N_);
  k2_=VectorAlloc(N_);
  k3_=VectorAlloc(N_);
  k4_=VectorAlloc(N_);
  k5_=VectorAlloc(N_);
  k6_=VectorAlloc(N_);
  temp1_=VectorAlloc(N_);

  x_[ns_]=time_;
  f_(x_, k1_);
  for(i=0; i<ns_; i++) {
    k1_[i]*=htime_;
    temp1_[i]=x_[i]+b21_*k1_[i];
  } /* for i */
  temp1_[ns_]=time_+a2_*htime_;
  f_(temp1_, k2_);
  for(i=0; i<ns_; i++) {
    k2_[i]*=htime_;
    temp1_[i]=x_[i]+b31_*k1_[i]+b32_*k2_[i];
  } /* for i */
  temp1_[ns_]=time_+a3_*htime_;
  f_(temp1_, k3_);
  for(i=0; i<ns_; i++) {
    k3_[i]*=htime_;
    temp1_[i]=x_[i]+b41_*k1_[i]+b42_*k2_[i]+b43_*k3_[i];
  } /* for i */
  temp1_[ns_]=time_+a4_*htime_;
  f_(temp1_, k4_);
  for(i=0; i<ns_; i++) {
    k4_[i]*=htime_;
    temp1_[i]=x_[i]+b51_*k1_[i]+b52_*k2_[i]+b53_*k3_[i]+b54_*k4_[i];
  } /* for i */
  temp1_[ns_]=time_+a5_*htime_;
  f_(temp1_, k5_);
  for(i=0; i<ns_; i++) {
    k5_[i]*=htime_;
    temp1_[i]=x_[i]+b61_*k1_[i]+b62_*k2_[i]+b63_*k3_[i]
      +b64_*k4_[i]+b65_*k5_[i];
  } /* for i */
  temp1_[ns_]=time_+a6_*htime_;
  f_(temp1_, k6_);
  for(i=0; i<ns_; i++) {
    k6_[i]*=htime_;
    xnew_[i]=x_[i]+c1_*k1_[i]+c3_*k3_[i]+c4_*k4_[i]+c6_*k6_[i];
    xerr_[i]=dc1_*k1_[i]+dc3_*k3_[i]+dc4_*k4_[i]+dc5_*k5_[i]+dc6_*k6_[i];
  } /* for i */
} /* RKFNCTime */


/***************************************************************************
  Solver_ = 7. (see also solver_ 4): For stiff non-autonomous systems, 
  [1, 3, 7, 9].

  Semi-implicit fourth order Rosenbrock method. The constants are due to
  Kaps and Rentrop (version GRK4T ,[3]), but have been transformed by the
  the method outlined in [1, pp. 738-742], [7], and [9].
  Adaptive stepsize control: Estimate of local error due to embedded third
  order method.
  The output from GRK4T: xnew_ and xerr_ contain respectively the new values
  and the estimate of the local error of the dependent variables.
****************************************************************************/

void GRK4TTime(int N_, int ns_, double time_, double htime_, double *x_, 
	       double *xnew_, double *xerr_, void (*f_)(double *, double *), 
	       void (*jac_)(double *)) {

  static const double gamma_    =   0.23100000000;
  static const double c21_      =  -5.07167544877;
  static const double c31_      =   6.02015272865;
  static const double c32_      =   0.159750684673;
  static const double c41_      =  -1.856343618677;
  static const double c42_      =  -8.50538085819;
  static const double c43_      =  -2.08407513602;
  static const double a21_      =   2.00000000000;
  static const double a31_      =   4.52470820736;
  static const double a32_      =   4.16352878860;
  /* a41_ = a31_, a42_ = a32_, a43_ = 0.0 */
  static const double m1_       =   3.95750374663;
  static const double m2_       =   4.62489238836;
  static const double m3_       =   0.617477263873;
  static const double m4_       =   1.282612945268;
  static const double dm1_      =  -2.30215540292;
  static const double dm2_      =  -3.07363448539;
  static const double dm3_      =   0.873280801802;
  static const double dm4_      =   1.282612945268;

  static const double AT2_      =   0.462000000000;
  static const double AT3_      =   0.880208333333;
  static const double BT1_      =   0.231000000000;
  static const double BT2_      =  -0.039629667552;
  static const double BT3_      =   0.550778939579;
  static const double BT4_      =  -0.055350984570;
 
  double **A_;
  int    *index_;
  double *k1_, *k2_, *k3_, *k4_;
  double *temp1_, *temp2_, *jactime_, dummy;
  int    i, j;

  A_=MatrixAlloc(N_);
  index_=IntVectorAlloc(N_);
  k1_=VectorAlloc(N_);
  k2_=VectorAlloc(N_);
  k3_=VectorAlloc(N_);
  k4_=VectorAlloc(N_);
  temp1_=VectorAlloc(N_);
  temp2_=VectorAlloc(N_);
  jactime_=VectorAlloc(N_);

  x_[ns_]=time_;
  dummy=1.0/(htime_*gamma_);
  jac_(x_);
  for(i=0; i<ns_; i++)
    jactime_[i]=htime_*jacobi_matx[i][ns_];
  for(i=0; i<ns_; i++)
    for(j=0; j<ns_; j++)
      A_[i][j]= -jacobi_matx[i][j];
  for(i=0; i<ns_; i++)
    A_[i][i]+=dummy;
  LUfact(ns_, A_, index_);

  f_(x_, k1_);
  for(i=0; i<ns_; i++)
    k1_[i]+=(BT1_*jactime_[i]);  
  LUsubst(ns_, A_, index_, k1_);

  temp1_[ns_]=time_+AT2_*htime_;
  for(i=0; i<ns_; i++)
    temp1_[i]=x_[i]+a21_*k1_[i];
  f_(temp1_, k2_);
  for(i=0; i<ns_; i++)
    k2_[i]+=((c21_*k1_[i])/htime_+BT2_*jactime_[i]);
  LUsubst(ns_, A_, index_, k2_);

  temp1_[ns_]=time_+AT3_*htime_;
  for(i=0; i<ns_; i++)
    temp1_[i]=x_[i]+a31_*k1_[i]+a32_*k2_[i];
  f_(temp1_, k3_);
  for(i=0; i<ns_; i++)  {
    k4_[i] = k3_[i];
    k3_[i]+=((c31_*k1_[i]+c32_*k2_[i])/htime_+BT3_*jactime_[i]);
  } /* for i */
  LUsubst(ns_, A_, index_, k3_);

  for(i=0; i<ns_; i++)
    k4_[i]+=((c41_*k1_[i]+c42_*k2_[i]+c43_*k3_[i])/htime_+BT4_*jactime_[i]);
  LUsubst(ns_, A_, index_, k4_);

  for(i=0; i<ns_; i++) {
    xnew_[i]= x_[i]+m1_*k1_[i]+m2_*k2_[i]+m3_*k3_[i]+m4_*k4_[i];
    xerr_[i]= dm1_*k1_[i]+dm2_*k2_[i]+dm3_*k3_[i]+dm4_*k4_[i];
  } /* for i */
  
  MatrixFree(N_, A_);
  VectorFree(N_, k1_);
  VectorFree(N_, k2_);
  VectorFree(N_, k3_);
  VectorFree(N_, k4_);
  IntVectorFree(N_, index_);
  VectorFree(N_, temp1_);
  VectorFree(N_, temp2_);
  VectorFree(N_, jactime_);
} /* GRK4TTime */
