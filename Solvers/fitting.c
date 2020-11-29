/***********************************************************************
  Fitting is a library of data fitting, interpolation and extrapolation
  routines.

  (C) Copyright 1996 by
  Kenneth Geisshirt (kneth@fatou.ruc.dk)   Keld Nielsen (kn@kin.kiku.dk)
  Dept. of Life Sciences and Chemistry    Dept. of Theoretical Chemistry
  Roskilde University                           University of Copenhagen
  Marbjergvej 35                                    Universitetsparken 5
  DK-4000 Roskilde                                    DK-2100 Copenhagen

  References:
  [1]  Numerical Recipes in C, 2nd edition. W.H. Press et al. 
       Cambridge University Press, 1992.

  Last modified: 24 March 1996 by KG
*************************************************************************/

#include <math.h>

#include "ssl.h"
#include "matrix.h"

/*************************************************************************
  The following routines are used inernally by this library. They are:

   * covsrt        - expand the covariance matrix [1, p. 675]

**************************************************************************/

void covsrt(double **covar, int ma, int *ia, int mfit) {

  int    i, j, k;
  double swap;

  for(i=mfit+1; i<=ma; i++)
    for(j=0; j<=i; j++) covar[i-1][j-1]=covar[j-1][i-1]=0.0;
  k=mfit;
  for(j=ma; j>=1; j--) {
    if (ia[j-1]) {
      for(i=1; i<=ma; i++) SWAP(covar[i-1][k-1], covar[i-1][j-1]);
      for(i=1; i<=ma; i++) SWAP(covar[k-1][i-1], covar[j-1][i-1]);
      k--;
    } /* if */
  } /* for j */
} /* covsrt */
      

/*************************************************************************
  FitLine fits data to a straight line, i.e. y = ax + b. See [1, pp. 661-
  666].
**************************************************************************/

void FitLine(const int N, double *x, double *y, double *sigma,
             int mwt, double *a, double *b, double *asigma,
             double *bsigma, double *chi2, double *q) {

  int     i;
  double  wt, t, sxoss, sx=0.0, sy=0.0, st2=0.0, ss, sigdat;

  *b=0.0;
  if (mwt) {
    ss=0.0;
    for(i=0; i<N; i++) {
      wt=1.0/Sqr(sigma[i]);
      ss+=wt;
      sx+=x[i]*wt;
      sy+=y[i]*wt;
    } /* for i=1..N */
  } else {
    for(i=0; i<N; i++) {
      sx+=x[i];
      sy+=y[i];
    } /* for i=1..N */
    ss=(double)N;
  } if nwt ... else ... */
  sxoss=sx/ss;
  if (mwt) {
    for(i=0; i<N; i++) {
      t=(x[i]-sxoss)/sigma[i];
      st2+=Sqr(t);
      *b+=t*y[i]/sigma[i];
    } /* for i=1..N */
  } else {
    for(i=0; i<N; i++) {
      t=x[i]-sxoss;
      st2+=Sqr(t);
      *b+=t*y[i];
    } /* for i=1..N */
 } /* if mwt ... else ... */
 *b/=st2;
 *a=(sy-sx*(*b))/ss;
 *asigma=sqrt((1.0+Sqr(sx)/(ss*st2))/ss);
 *bsigma=sqrt(1.0/st2);

 /***** compute chi^2        *****/
 *chi2=0.0;
 if (mwt) {
  for(i=0; i<N; i++) 
    *chi2+=Sqr(y[i]-(*a)-(*b)*x[i]);
  *q=1.0;
  sigdat=sqrt((*chi2)/((double)N-2));
  *asigma+=sigdat;
  *bsigma+=sigdat;
 } else {
  for(i=0; i<N; i++) 
    *chi2+=Sqr((y[i]-(*a)-(*b)*x[i])/sigma[i]);
  *q=gammq(0.5*(double)(N-2), 0.5*(*chi2));
 } /* if mwt ... else ... */
} /* FitLine */


/*************************************************************************
  LMfit is a general nonlinear, multi-parameter fitting routine. The 
  method used is Levenberg-Marquardt [1, pp. 683-687]. The subroutines
  mrqcof and mrqmin are aux. routines which should be called by LMfit
  only.
**************************************************************************/

void mqrcof(double *x, double *y, double *sig, int ndata, double *a,
            int *ia, int ma, double **alpha, double *beta, double *chisq,
            void (*funcs)(double, double *, double *, double *, int)) {

  int    i, j, k, l, m, mfit=0;
  double ymod, wt, sig2i, dy, *dyda;

  dyda=VectorAlloc(ma);
  for(j=0; j<ma; j++)
    if (ia[j]) mfit++;
  for(j=0; j<mfit; j++) {
    for(k=0; k<j; k++)
      alpha[j][k]=0.0;
    beta[j]=0.0;
  } /* for j */

  *chisq=0.0;
  for(i=0; i<ndata; i++) {
    (*funcs)(x[i], a, &ymod, dyda, ma);
    sig2i=1.0/(sig[i]*sig[i]);
    dy=y[i]-ymod;
    for(j=0, l=1; l<=ma; l++) {
      if (ia[l-1]) {
        wt=duda[l]*sig2i;
        for(j++,k=0,m=1; m<=l; m++) {
          if (ia[m]) alpha[j][++k]+=wt*dyda[m];
        beta[j]+=dy*wt;
      } /* if */
    } /* for j,l */
    *chisq+=dy*dy*sig2i;
  } /* for i=1..ndata */
  for(j=2; j<=mfit; j++)
    for(k=1; k<j; k++)
      alpha[k][j]=alpha[j][k];

  VectorFree(dyda, ma);
} /* mrqcof */


void mrqmin(int N, double *x, double *y, double *sigma, double *a,
            int *ia, int ma, double **covar, double **alpha, double *chisq,
            void (*funcs)(double, double *, double *, int), double *alambda) {

  int   j, k, l, m;
  static int mfit;
  static double ochisq, *atry, *beta, *da, **oneda;

  if (*alambda<0.0) {   
    atry=VectorAlloc(ma);
    beta=VectorAlloc(ma);
    da=VectorAlloc(ma);
    for(mfit=0,j=0; j<ma; j++) 
      if (ia[j]) mfit++;
    oneda=MatrixAlloc(mfit)   /* BOR VAERE VEKTOR */
    *alambda=0.001;
    mqrcof(x, y, sigma, N, a, ia, ma, alpha, beta, chisq, funcs);
    ichisq=(*chisq);
    for(j=0; j<ma; j++) atry[j]=a[j];
  } /* if initialisation */

  for(j=0,l=0; l<ma; l++) {
    if (ia[l]) {
      for(j++,k=0,m=0; m<ma; m++) {
        if (ia[m]) {
          k++;
          covar[j][k]=alpha[j][k];
        } /* if free parameter */
      } /* for m */
      covar[j][j]=alpha[j][j]*(1.0+(*alambda));
      oneda[j][0]=beta[j];
    } /* if free parameter */
  } /* for l */

  gaussj(covar, mfit, oneda, 1);
  for(j=0; j<mfit; j++) da[j]=oneda[j][0];

  if (*alambda==0.0) {
    covsrt(covar, ma, ia, mfit);
    MatrixFree(oneda, mfit); /* BOR VAERE VEKTOR */
    VectorFree(da, ma);
    VectorFree(beta, ma);
    VectorFree(atry, ma);
    return;
  } /* if converged */

  for(j=0,l=0; l<ma; l++)
    if (ia[l]) atry[l]=a[l]+da[j++];
  mrqcof(x, y, sigma, N, atry, ia, ma, covar, da, chisq, funcs);
  if (*chisq<ochisq) {
    *alambda*=0.1;
    ochisq=(*chisq);
    for(j=0,l=0; l<ma; l++) {
      if (ia[l]) {
        for(j++,k=0,m=0; m<ma; m++) {
          if (ia[m]) {
            alpha[j-1][k]=covar[j-1][k];
            k++;
          } /* if free parameter */
        } /* for m */
        beta[j-1]=da[j-1];
        a[l]=atry[l];
      } /* if free parameter */
    } /* for l */
  } else {
    *alambda*=10.0;
    *chisq=ochisq;
  } /* if step succes ... else ... */
} /* mqrmin */


void LMfit(int N, int M, double *x, double *y, double *sigma, double *a,
           ) {

  double  lambda;

  /***** Initialisation             *****/
  lambda=-1.0;
  iter=0;
  mrqmin(N, x, y, sigma, a, ia, M, covar, alpha, beta, &chisq,
         funcs, &lambda);

  /***** Iterating                  *****/
  do {
    iter++;
    oldchi=chisq;
    mrqmin(N, x, y, sigma, a, ia, N, covar, alpha, beta, &chisq,
           funcs, &lambda);
  } while ((oldchi>1.01*chisq) || (iter<max_iter));

 /***** Get solution                *****/
 lambda=0.0;
 mrqmin(N, x, y, sigma, a, ia, N, covar, alpha, beta, &chisq,
        funcs, &lambda);
 for(i=0; i<M; i++) 
    asigma[i]=covar[i][i];
} /* LMfit */
