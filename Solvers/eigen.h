/*****************************************************************************
  Eigen is a library for computing eigenvalues and eigenvectors of general
  matrices. The main routine is Eigen; ComplexVectorAlloc, ComplexMatrixAlloc
  and the matching Free functions manage the arrays it fills in.

  The meaning of the arguments to Eigen is:
    1.   The dimension of the general matrix (n).
    2.   A general matrix (A).   
    3.   The maximal number of iterations.
    4.   The precision.
    5.   A vector with the eigenvalues.
    6.   A matrix with the eigenvectors.


  CopyWrong 1994 by
  Kenneth Geisshirt (kneth@osc.kiku.dk)
  Department of Theoretical Chemistry
  H.C. Orsted Institute
  Universitetsparken 5
  2100 Copenhagen
  Denmark

  Last updated: 7 February 1995 , by KN
*****************************************************************************/

#ifndef _EIGEN_LIB_
#define _EIGEN_LIB_

#include <stdio.h>
#include <complex.h>

extern void Eigen(int, int, double **, int, double, int, double complex *,
                  double complex **);
extern double complex  *ComplexVectorAlloc(const int);
extern double complex **ComplexMatrixAlloc(const int);
extern void ComplexVectorFree(const int, double complex *);
extern void ComplexMatrixFree(const int, double complex **);

#endif
