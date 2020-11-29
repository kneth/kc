/*************************************************************************
  (C) Kenneth Geisshirt (kneth@fatou.ruc.dk)
      Department of Life Sciences and Chemistry
      Roskilde University
      P.O. Box 260
      4000 Roskilde
      Denmark

   This is header file of a small library for sparse matrices.

   Last updated: 18 April 1995
**************************************************************************/

#ifndef _SPARSE_MATRIX_LIB_
#define _SPARSE_MATRIX_LIB_

struct sparse_matrix {
  double  *A;    /* the elements */
  int     *CNR;
  int     *RNR;
  int     nz;    /* number of elements */
  int     N;     /* number of rows/columns */


