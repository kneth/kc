/***************************************************************************
  The source file contains the driver for solving ODEs.

  CopyWrong 1994-1995 by
  Kenneth Geisshirt (kneth@fatou.ruc.dk) 
  Department of Life Sciences and Chemistry
  Roskilde University
  P.O. Box 260
  4000 Roskilde
  Denmark

  and 
  Keld Nielsen (kn@kiku.dk)
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

  Last updated: 5 June 1995 by KN
****************************************************************************/

#define VERSION_ "1.05"

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <signal.h>

/***************************************************************************
  The solvers for ODEs and the service routines are imported below.
****************************************************************************/

#include "odesolv.h"
#include "odeserv.h"


/***************************************************************************
  The matrix manipulating routines like linear equation solvers are all
  implemented in a small library. It is imported below. 
  
  Routines for integrating functions are also imported below.
****************************************************************************/

#include "matrix.h"

/***************************************************************************
  The program has a number of global variables. Each variable has an 
  underscore as the last character in order to eliminate conflicts with 
  names in the model. 

  time_       The independent variable of the ODEs 
  stime_      Initial value for the independent variable. Default: 0 
  etime_      Final value for the independent variable. Default: 200 
  dtime_      Output at equidistant points. Also functioning as the largest 
              possible stepsize. Default: 2.0 
  htime_      Stepsize 
  epsr_       Relative error tolerance. Default: 1.0E-5
  epsa_       Absolute error tolerance. Default: 1.0E-15 
  epsmin_     Relative machine precision. Default: 1.0E-16 
  timenew_    New value of the independent variable before a new step of 
              integration is accepted. 
  htimenew_   New value of the stepsize before a new step of integration is 
              accepted. 
  errlim_     The ratio between the estimate of the (scaled) local truncation 
              error and the relative tolenrance. 
	      If errlim_>1.0 the steplength (hstep_) is rejected otherwise 
	      the step of integration is accepted. 
  thtime_     The value of the independent variable for the next output.
  step_adjustment_factor 
              Safety factor in the calculation of a new stepsize. Default: 0.9 
  order_of_method
              The order of the integration method 
  steplimit_increase
              Upper bound for the ratio of two consecutive steps. Default: 1.5
  steplimit_decrease
              Lower bound for the ratio of two consecutive steps. Default: 0.5
  step_increase_factor == -1/order_of_method. 
	      Used in the calculation of the new steplength. 
  errcon_     Smallest value of errlim_, [1] 
  htimemin_   The minimal allowed stepsize. Default: Relative machine 
              precision.

  mainmode_   1==ODEs, 2==PDEs, 3==compute LCEs. Default: 1
  solver_     Method of integration. Default: 1 
  prnmode_    Mode of output. 0: Equidistant and extrema points. 
                              1: Only equidistant points. 
                              2: Quantum chemistry mode (don't use).
              Default: 0 
  scaling_    Method of scaling. 0: Strict relative error [8], 
                                 1: Scaling device in error estimate due 
				 to Deuflhard et al., [8].
  debug_       Mode of debug: 0: No output during numerical integration.
                              1: Time, steplength and values of the variables
			         are printed.
			      2: Different control values are printed. 
			         Extension of mode 1.
			      3: Initial values of different control values are
			         printed. 
               Default: 0
  finish_      If the user interupt the program, finish_ will be set to
               1, and the program will be terminated nicely. See also the 
	       signal handler.
  pert_        Keeps track of pertubation of the differential equations.
  ptime_       Time for a perturbation.
  dptime_      Time between perturbations.
  no_pert      Number of perturbations done so far.
  next_pert    Time for next perturbation.
  outfile      Handler to output file.
  i            Simple counter.
  datafilename_Name of output file.

  equa         Number of differential equations.
  ns_ = equa-1 Number of dependent variables in a nonautonomous 
               differential equations. 
  xnew_        New values of the dependent variables.
  xerr_        Estimate of the unscaled local truncation errors for the 
               dependent variables. 
  begin_print  Time for first print out.
  species_     A table with the names of the independent variables.


  The file model.c contains the C-code that depends on a specific model. 
  It is generated by kc, see [4].
*****************************************************************************/

static double time_, stime_, etime_, dtime_, htime_, epsr_, epsa_;
static double epsmin_, timenew_, htimenew_, errlim_, thtime_, ptime_;
static double step_adjustment_factor, order_of_method, steplimit_increase;
static double steplimit_decrease, step_increase_factor, errcon_, htimemin_;
static int    mainmode_, solver_, prnmode_, scaling_, pert_=1, i;
static double dptime_, next_pert, no_pert;
static double begin_print;
static FILE   *outfile;
static char   datafilename_[35];
static int    iter_, debug_=1, finish_=0;
static int    steprejection_=0;         /* Number of step rejections */
static double cmax_, cmin_;
int    nospecieserr_;    /* Number of variables used in the estimate
			    of local error */

#include "model.c" 

static const int       ns_ = equa-1;
static double xnew_[equa], xerr_[equa];
static double fx_[equa], temp_[equa], xdt_[equa];

/****************************************************************************
  The routine UserInterrupt catch the interrupt signal from the operating
  system. This signal is sent when the user wants to abort the program
  before normal termination. The feature ensure that the buffers are flushed,
  and the termination in general is done "nicely". The signal handler is 
  very simple - the actual work is done by the main program.
*****************************************************************************/

void UserInterrupt(int dummy) {

  finish_=1;
} /* UserInterrupt */


/***************************************************************************
  PrintDebug is the debugging routine of the program. It prints the 
  information specified by the variable debug_. The information is printed
  on the standard output.
****************************************************************************/

void PrintDebug(int flag_, double t1_, double ht1_, double *x1_, double t2_,
		double ht2_, double *x2_) {

  int i;

  switch (flag_) {
  case 0:
    /* do nothing */
    break;
  case 1:
    printf("\33HTime: %e        Steplength: %e\33K\n\33K\n", t2_, ht2_);
    for(i=0; i<equa; i++)
      printf("%-12s%e\33K\n", species[i], x2_[i]);
    printf("\33J");
    break;
  case 2:
    printf("Accept point if errlim<=1.0: errlim=%e\t(errcon=%e)\n\n",
	   errlim_, errcon_);
    printf("Time: %e\tNew time: %e\t(Next print out: %e)\n", t1_, t2_,
	   thtime_);
    printf("Steplth: %e\tNew steplth: %e\t\n\n", ht1_, ht2_);
    for(i=0; i<equa; i++)
      printf("%s\t%e\t\t%e\t(xerr, xscal: %e  %e)\n\n", species[i], x1_[i],
	     x2_[i], fabs(xerr_[i]), xscal_[i]);
    printf("Decrease of rejected step, if steprejection=2: %d\n\n\n\n", 
	   steprejection_);
    break;
  case 3:
    printf("****************** Initial Values *********************\n\n");
    printf("epsr: %e\tepsa: %e\tepsmin: %e\terrcon: %e\n\n", epsr_, epsa_,
	   epsmin_, errcon_);
    printf("time: %e\tdtime: %e\tetime: %e\n\n", time_, dtime_, etime_);
    printf("steplength: %e\thtimemin: %e\n\n", htime_, htimemin_);
    printf("stepadjust: %e\tmaxinc: %e\tmininc: %e\n\n", 
	   step_adjustment_factor, steplimit_increase, steplimit_decrease);
    printf("method: %d\torder of method: %3.0f\n\n", solver_, order_of_method);
    printf("No. of differentail equations: %d\n\n", equa);
    printf("Initial values of the variables:\n");
    for(i=0; i<equa; i++)
      printf("  %s\t%e\n", species[i], x_[i]);
    printf("\nscaling method: %d\tcmax: %e\tcmin: %e\n\n", scaling_, cmax_, 
	   cmin_);
    printf("Initial scaling values:\n");
    for(i=0; i<equa; i++)
      printf("  %s\t%e\n", species[i], xscal_[i]);
    printf("\nmainmode: %d\tprnmode: %d\tdebug: %d\tname of datefile: %s\n\n",
	   mainmode_, prnmode_, debug_, datafilename_);
    printf("****************** Initial Values *********************\n\n");
    break;
  } /* switch flag */
} /* PrintDebug */


int ChangeSignVectors(int N_, double *x_, double *xnew_, double *xmax_,
		      double *xmin_, int *osg_, int *nsg_) {

  int        i, flag_=0;
  static int first_time=1;

  for(i=0; i<N_; i++) {
    if ((xnew_[i]-x_[i])<0.0)
      nsg_[i]=-1;
    else if ((xnew_[i]-x_[i])>0.0)
      nsg_[i]=1;
    else 
      nsg_[i]=0;
  } /* for i */

  if (first_time) {
    first_time=0;
    for(i=0; i<N_; i++)
      osg_[i]=nsg_[i];
    return flag_;
  } /* if */
  
  for(i=0; i<N_; i++) {
    if ((osg_[i]!=nsg_[i]) && (nsg_[i]==-1)) {
      osg_[i]=nsg_[i];
      if (fabs(x_[i]-xmin_[i])>epsa_)
	flag_++;
    } /* if */
    if ((osg_[i]!=nsg_[i]) && (nsg_[i]==1)) {
      osg_[i]=nsg_[i];
      if (fabs(x_[i]-xmax_[i])>epsa_)
	flag_++;
    } /* if */
  } /* for i */
  return flag_;
} /* ChangeSignVectors */


/***************************************************************************
  Main initialises the system, and it is also the driver routine for the
  numerical schemes implemented.

  The idea is that each solver takes one time step, and main contains a 
  loop going from time=initial time to termination time. The step length
  controller is implemented as part of main.
****************************************************************************/

void main(void) {
  
  double xmin_[equa], xmax_[equa], hlp_;
  time_t timer_;           /* time and date */
  int    i, j, dyn_, csv_, prndt_;
  int    nsg_[equa], osg_[equa];
  

  printf("KKsolver v%s, CopyWrong 1994-1995 by Keld Nielsen and Kenneth Geisshirt\n", VERSION_);

  /* Find time and date */
  timer_=time(&timer_);

  /* Set up signal handler */
  (void) signal(SIGINT, &UserInterrupt);

  /* Allocate the jacobian matrix */
  jacobi_matx=MatrixAlloc(equa);

  /* Initialization of the model dependent variables. */
  InitValues();

  /* Initialize parameters for various printing modes */
  switch (prnmode_) {
  case 0:
    for(i=0; i<equa; i++) {
      nsg_[i]=osg_[i]=0;
      xmin_[i]=xmax_[i]=0.0;
    } /* for i */
    prndt_=0;
    break;
  } /* switch prnmode */

  /* Setup "precision" */
  if (epsmin_==0.0)
    epsmin_=FindMachEps();

  /* Time equal to the initial time. */
  time_=stime_;
  
  /* The next time a output will be performed */
  thtime_=time_+dtime_;

  /* Program aborted: If initial time is larger than the end time */
  if (time_>etime_)
    kcerror("stime > etime");

  /* Program aborted: If initial time is less than the zero */
  if (time_<0.0)
    kcerror("stime < 0.0");

  /* Program aborted: If initial stepsize is less than the zero */
  if (htime_<0.0)
    kcerror("htime < 0.0");

  /* Program aborted: If requested output at eqvidistant points 
     is not resonable */
  if (dtime_<0.0)
    kcerror("dtime < 0.0");

  /* Smallest absoulte error not less than relative machine precision */
  if (epsa_<epsmin_)
    epsa_=epsmin_;

  /* Smallest stepsize not less than relative machine precision */
  if (htimemin_<epsmin_)
    htimemin_=epsmin_;

  /* Smallest initial stepsize equal to htimemin_ */
  if (htime_<htimemin_)
    htime_=htimemin_;

  /* Largest initial stepsize equal to 0.5*dtime_  */
  if (htime_>=0.5*dtime_)
    htime_=0.5*dtime_;

  /* Opening output file */
  outfile=fopen(datafilename_, "w");
  fprintf(outfile, "# Output from kksolver v%s, %s", VERSION_, 
	  ctime(&timer_));
  fprintf(outfile, "# CopyWrong 1994-1995 by Keld Nielsen and Kenneth Geisshirt\n");
  dyn_=1;
  for(i=0; i<equa; i++) 
    if (do_print[i]==1) {
      fprintf(outfile, "# Dynamical var. no. %d: %s\n", dyn_, species[i]);
      dyn_++;
    } /* if */
  
#ifdef _DO_PERT_
  no_pert=0.0;
  next_pert=ptime_;
#endif

  /* Scaling method: Initialization */
  if (scaling_==0) {               /* Strict relative error */
    for(i=0; i<equa; i++)
      xscal_[i]=MaxPair(epsa_, fabs(x_[i]));
  } else {                         /* Scaling due to Deuflhard et al. ([8]) */
    cmax_= fabs(x_[0]);
    for (i=0; i<equa; i++)
      cmax_=MaxPair(cmax_, fabs(x_[i]));
    cmin_= cmax_*epsmin_/epsr_;
    for (i=0; i<equa; i++)
       xscal_[i]=MaxPair(cmin_, fabs(x_[i]));
  } /* if else */

  /* Initialising the various schemes */
  switch (solver_) {
  case 1:                          /* Calahan - stiff, autonom        */
    order_of_method=3.0;
    nospecieserr_=equa;
    break;
  case 2:                          /* GRK4T - stiff, autonom          */
    order_of_method=4.0;
    nospecieserr_=equa;
    break;
  case 3:                          /* RKFNC - non stiff, autonom      */
    order_of_method=5.0;
    nospecieserr_=equa;
    break;
  case 4:                         /* RK4 - non stiff, autonom         */
    order_of_method=4.0;
    nospecieserr_=equa;
    break;
  case 5:                         /* GRK4TTime - stiff, non autonom   */
    order_of_method=4.0;
    nospecieserr_=equa-1;
    break;
  case 6:                         /* RKFNCTime-non stiff, non autonom */
    order_of_method=5.0;
    nospecieserr_=equa-1;
    break;
  } /* switch */

  step_increase_factor=(-1.0/(order_of_method));

  /* Smallest value of errlim_, see [1] */
  errcon_=pow(steplimit_increase/step_adjustment_factor,
	      -1.0*(order_of_method+1.0));

  /* Print out of initial values */
  if (time_>=begin_print)
    PrintState(nospecieserr_, prnmode_, 0, time_, x_, outfile, do_print);

  if (debug_==3)
    PrintDebug(debug_, time_, htime_, x_, time_, htime_, x_);

  /* Main part of the integration algorithm */
  while ((time_<etime_) && (finish_==0)) {

    /* Integration of one step */
    switch (solver_) {
    case 1:                               /* Calahan       */
      Calahan(equa, htime_, x_, xnew_, xerr_, &reac, &jacobi);
      break;
    case 2:                               /* GRK4T         */
      GRK4T(equa, htime_, x_, xnew_, xerr_, &reac, &jacobi);
      break;
    case 3:                               /* RKFNC         */
      RKFNC(equa, htime_, x_, xnew_, xerr_, &reac);
      break;
    case 4:                               /* RK4           */
      RK4(equa, htime_, x_, xnew_, xerr_, &reac);
      break;
    case 5:                               /* GRK4Time      */
      GRK4TTime(equa, ns_, time_, htime_, x_, xnew_, xerr_, &reac, jacobi);
      break;
    case 6:                               /* RKFNCTime     */
      RKFNCTime(equa, ns_, time_, htime_, x_, xnew_, xerr_, &reac);
      break;
    } /* switch (solver_) */

    if (scaling_==0) {
      for(i=0; i<equa; i++)
	xscal_[i]=MaxPair(epsa_, fabs(xnew_[i]));
    } else {
      for(i=0; i<equa; i++)
	xscal_[i]=MaxPair(xscal_[i], MaxPair(cmin_, fabs(x_[i])));
    } /* else */

    /* Scaling the error estimate */
    for(i=0; i<equa; i++)
      temp_[i]=xerr_[i]/xscal_[i];

    /* Maximal values of (local error estimate/relative error tolerance) */
    errlim_=MaxVec(nospecieserr_, temp_)/epsr_;

    /* New stepsize: The new stepsize is determined from: SAF*(EST/RTOL)**(SIF)
		     (SAF=step_adjustment_factor, EST=local error estimate,
		     RTOL=relative error tolerance, **(SIF)=pow of 
		     step_increase_factor).
		     The new stepsize is restricted to the following values:
		       [steplimit_decrease*h_old;steplimit_increase*h_old].
		     If the new stepsize is larger than dtime_, then it 
		     is changed to dtime_. */
    
    if (errlim_<errcon_) {
      htimenew_=MinPair(dtime_, htime_*steplimit_increase); /* Small values of
							    errlim_, see [1] */
    } else {
      htimenew_ = MinPair(dtime_, htime_*
			  MaxPair(steplimit_decrease,
				  MinPair(steplimit_increase, 
				  step_adjustment_factor
				  *pow(errlim_, step_increase_factor))));
    } /* else */


   /* Program aborted: If stepsize less than the smallest 
                       acceptable stepsize */
    if (htimenew_<htimemin_)
      kcerror("htimenew < htimemin");

    /* Integration does not exceed etime_ */
    if (thtime_>etime_)
      thtime_=etime_;

    if (debug_==2)
      PrintDebug(debug_, time_, htime_, x_, time_+htime_, htimenew_, xnew_);

    if (errlim_>1.0) {  /* The new stepsize is not accepted if the 
			   local error estimate is larger than the 
			   relative error tolerance */
      htime_=htimenew_;
      if (steprejection_!=2) {
	steprejection_++;
      } else {          /* If the new stepsize has been rejected more 
			   than two times the new stepsize is drastical 
			   reduced; h_new= h_old/10.
			   Method suggested by Hairer and Wanner [9] */
	steprejection_=0;
	htime_= htimenew_/10.0;
      } /* else */
    } else  {          /* The new stepsize is accepted if the local error 
			  estimate is less or equal to the relative error 
			  tolerance */
      timenew_=time_+htime_;
      steprejection_=0;

    if (timenew_<ptime_) {
      if (timenew_<thtime_) {  /* The new time do not exceed the time
				  of the next output */
	if (prnmode_==0) {     /* Output: The extrema values */
	  csv_=ChangeSignVectors(nospecieserr_, x_, xnew_, xmin_, xmax_,
				 osg_, nsg_);
	  if ((csv_) && (prndt_==0) && (time_>=begin_print)) 
	    PrintState(nospecieserr_, prnmode_, 2, time_, x_, outfile, 
		       do_print);
	  if (prndt_)
	    prndt_=0;
	} /* if */
	time_=timenew_;
	htime_=htimenew_;
	for(i=0; i<equa; i++)
	  x_[i]=xnew_[i];
      } else {
	if (timenew_==thtime_) { /* Output: the new time is equal to the
				    time for the next output */
	  if (prnmode_==0) {
	    csv_=ChangeSignVectors(nospecieserr_, x_, xnew_, xmax_, xmin_, 
				   osg_, nsg_);
	    prndt_=1;
	  } /* if */
	  time_=timenew_;
	  htime_=htimenew_;
	  for(i=0; i<equa; i++)
	    x_[i]=xnew_[i];
	  thtime_+=dtime_;
	  if (time_>=begin_print)
	    PrintState(nospecieserr_, prnmode_, 1, time_, x_, outfile, 
		       do_print);
	  if (debug_==1)
	    PrintDebug(debug_, timenew_, htimenew_, xnew_, time_, htime_, x_);
	} else  {              /* The new time exceed the time for the 
				  next output. The differential equations 
				  are integrated from (t) to (t+dt), and 
				  the solution at t+dt is printed. The 
				  integration begin again from (t+h, xnew_) */
	  hlp_=thtime_-time_;
	  switch(solver_) {
	  case 1:                               /* Calahan       */
	    CalahanOneStep(equa, hlp_, x_, xdt_, &reac, &jacobi);
	    break;
	  case 2:                               /* GRK4T         */
	    GRK4T(equa, hlp_, x_, xdt_, xerr_, &reac, &jacobi);
	    break;
	  case 3:                               /* RKFNC         */
	    RKFNC(equa, hlp_, x_, xdt_, xerr_, &reac);
	    break;
	  case 4:                               /* RK4           */
	    RK4OneStep(equa, hlp_, x_, xdt_, &reac);
	    break;
	  case 5:                               /* GRK4TTime     */
	    GRK4TTime(equa, ns_, time_, hlp_, x_, xdt_, xerr_, &reac, &jacobi);
	    break;
	  case 6:                               /* RKFNCTime     */
	    RKFNCTime(equa, ns_, time_, hlp_, x_, xdt_, xerr_, &reac);
	    break;
	  } /* switch (solver_) */

	  if (prnmode_==0) {     /* Output: The extrema values */
	    csv_=ChangeSignVectors(nospecieserr_, x_, xnew_, xmax_, xmin_,
				   osg_, nsg_);
	    prndt_=1;
	  }
	  if (debug_==1)
	    PrintDebug(debug_, time_, htime_, x_, thtime_, htime_, xdt_);
	  time_=thtime_;
	  htime_=htimenew_;
	  thtime_+=dtime_;
	  for(i=0; i<equa; i++)
	    x_[i]=xdt_[i];
	  if (time_>=begin_print)
	    PrintState(nospecieserr_,prnmode_,1,time_,x_,outfile,do_print);
	} /* else */
      } /* else */
  } /* ptime_>timenew_ */ 
  else {
    if (timenew_>ptime_) {
  PrintState(nospecieserr_, prnmode_, 1, time_, x_, outfile, do_print);
  hlp_=thtime_-time_;
  switch(solver_) {
  case 1:                               /* Calahan       */
    CalahanOneStep(equa, hlp_, x_, xdt_, &reac, &jacobi);
    break;
  case 2:                               /* GRK4T         */
    GRK4T(equa, hlp_, x_, xdt_, xerr_, &reac, &jacobi);
    break;
  case 3:                               /* RKFNC         */
    RKFNC(equa, hlp_, x_, xdt_, xerr_, &reac);
    break;
  case 4:                               /* RK4           */
    RK4OneStep(equa, hlp_, x_, xdt_, &reac);
    break;
  case 5:                               /* GRK4TTime     */
    GRK4TTime(equa, ns_, time_, hlp_, x_, xdt_, xerr_, &reac, &jacobi);
    break;
  case 6:                               /* RKFNCTime     */
    RKFNCTime(equa, ns_, time_, hlp_, x_, xdt_, xerr_, &reac);
    break;
  } /* switch (solver_) */
  time_=thtime_;
  htime_=htimenew_;
  thtime_+=dtime_;
  for(i=0; i<equa; i++)
    x_[i]=xdt_[i];
  fprintf(outfile, "# Perturbation at time=%e\n", time_);
  if (time_>=begin_print)
    PrintState(nospecieserr_,prnmode_,1,time_,x_,outfile,do_print);
  for(i=0; i<equa; i++)
    x_[i]+=x_pert[i];
  if (time_>=begin_print)
    PrintState(nospecieserr_,prnmode_,1,time_,x_,outfile,do_print);
  ptime_ += etime_+1.0;
    } else { 
    if (timenew_=ptime_) {
  time_=timenew_;
  htime_=htimenew_;
  for(i=0; i<equa; i++)
    x_[i]=xnew_[i];
  thtime_+=dtime_;
  fprintf(outfile, "# Perturbation at time=%e\n", time_);
  if (time_>=begin_print)
    PrintState(nospecieserr_, prnmode_, 1, time_, x_, outfile, do_print);
  for(i=0; i<equa; i++)
    x_[i]+=x_pert[i];
  if (time_>=begin_print)
    PrintState(nospecieserr_,prnmode_,1,time_,x_,outfile,do_print);
  ptime_ += etime_+1.0;
      }
    }
   }

    } /* errlim_ <= 1.0 */

  } /* while (time_<etime_) */

  if (finish_==1)
    fprintf(outfile, "# The user interrupted the program at time=%e\n", time_);
  fclose(outfile);

  MatrixFree(equa, jacobi_matx);
} /* main */
