/***************************************************************************
  This source file contains service routines for the solvers in odesolv.c

  CopyWrong 1994-1996 by
  Kenneth Geisshirt (kneth@fatou.ruc.dk)       Keld Nielsen (kn@kin.kiku.dk)
  Dept. of Life Sciences and Chemistry        Dept. of Theoretical Chemistry
  Roskilde University                               University of Copenhagen
  P.O. Box 260                                          Universitetsparken 5
  4000 Roskilde                                         2100 Copenhagen East
  Denmark                                                            Denmark

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
  [11] J.D. Smith. Design and Analysis of Algorithms. PWS-Kent, 1989.

  Last updated: 15 February 1995 by KG
****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <time.h>

/****************************************************************************
  The service routines are:

  o kcerror     - prints an error message and aborts the program.  
  o FindMachEps - find the machine's precision
  o MaxVec      - find the (abs) largest number in a vector
  o MaxPair     - find the largest of two numbers
  o MinPair     - find the smallest of two numbers
  o PrintState  - print the independent variable (time_) and the dependent
                  variables (x_). The data is written to the disk. The file
                  handler is outfile. The output is in ASCII characters, and
                  can be found in the following way:
                    time_  x_[0] x_[1]  ...  x_[equa-1] x_[equa]
                  Three modes of output. prnmode_=0: Equidistant and extrema
                  points printed, and prnmode_=1: Only equidistant points.
                  Default: 0.
  o BSort       - Bubble sort. The sorted list is a vector of real numbers,
                  [11, p. 157].
  o KronDelta   - Kronecker delta.
*****************************************************************************/

void kcerror(char *str) {

  fprintf(stderr, "Failure in numerical integration: %s.\n", str);
  fprintf(stderr, "Aborting the program - sorry.\n");
  exit(-1);
} /* kcerror */


double FindMachEps(void) {

  double temp=1.0;

  while ((1.0+temp)!=1.0)
    temp/=10.0;
  return temp;
} /* FindMachEps */


double MaxVec(int nosp_, double *x_) {

  double temp1_=0.0, temp2_;
  int    i_;

  for(i_=0; i_<nosp_; i_++) {
    temp2_=fabs(x_[i_]);
    if (temp2_>temp1_)
      temp1_=temp2_;
  } /* for i_ */
  return temp1_;
} /* MaxVec */


double MaxPair(double x_, double y_) {

  if (x_>y_)
    return x_;
  else
    return y_;
} /* MaxPair */


double MinPair(double x_, double y_) {

  if (x_<y_)
    return x_;
  else
    return y_;
} /* MinPair */


void PrintState(int N_, int prnmode_, int prnflag_, double time_, 
		double *x_, FILE *outfile_, int *do_prn) {
  
  int    i_, j_, flag_;
  static int first_time = 1;

  switch (prnmode_) {
  case 0:                            /* printing extrema */
    switch (prnflag_) {
    case 0:
      fprintf(outfile_, "# Initial values of the variables\n");
      break;
    case 1:
      break;
    case 2:
      fprintf(outfile_, "# Minimum or maximum of a variable\n");
      break;
    } /* switch prnflag_ */
    fprintf(outfile_, "%e ", time_);
    for(i_=0; i_<N_; i_++)
      if (do_prn[i_]==1)
	fprintf(outfile_, "%e ", x_[i_]);
    fprintf(outfile_, "\n");
    break;
  case 1:
    fprintf(outfile_, "%e ", time_);
    for(i_=0; i_<N_; i_++)
      if (do_prn[i_]==1)
	fprintf(outfile_, "%e ", x_[i_]);
    fprintf(outfile_, "\n");
    break;
  } /* switch prnmode */
} /* PrintState */

void BSort(const int N, double *list) {

  double   swap;
  int      i;

  for(i=0; i<(N-1); i++) 
    if (list[i]<list[i+1]) {
      swap=list[i];
      list[i]=list[i+1];
      list[i+1]=swap;
    } /* if */
} /* BSort */

double KronDelta(int i, int j) {

  if (i==j)
    return 1.0;
  else
    return 0.0;
} /* KronDelta */

