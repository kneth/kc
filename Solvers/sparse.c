/*************************************************************************
  (C) Kenneth Geisshirt (kneth@fatou.ruc.dk)
      Department of Life Sciences and Chemistry
      Roskilde University
      P.O. Box 260
      4000 Roskilde
      Denmark

   This is a simple implementation of sparse matrices.

   References:
   [1] Direct Methods for Sparse Matrices. O. Osterby and Z. Zlatev.
       Lecture Notes in Computer Science, 157, 1983. Springer-Verlag.
   [2] Computational Methods for General Sparse Matrices. Z. Zlatev,
       Kluwer Academic, 1991.
   [3] Sparse Matrix Techniques. Eds. A. Dold and B. Eckmann. Lecture
       Notes in Mathematics, 572, 1979. Springer-Verlag.

   Last updated: 19 April 1995
**************************************************************************/

#include <malloc.h>
#include <stdio.h>

#include "sparse.h"


/*************************************************************************
  SPorder reorders the data stored in the matrix, [1, pp. 20-21
**************************************************************************/

void SPorder(SPmatrix A) {
  
  int   i, j;
  
  for(i=0; i<A->N; i++) {
    A->ha[i][0]=0;
    A->ha[i][2]=0;
    A->ha[i][3]=0;
    A->ha[i][5]=0;
  } /* for i */
  for(i=0; i<A->NZ; i++) {
    j=A->cnr[i];
    A->ha[j][6]++;
    j=A->rnr[i];
    A->ha[j][3]++;
  } /* for i */
  n1=A->N-1;
  for(i=0; i<n1; i++) {
    A->ha[i+1][0]=A->ha[i][1]+A->ha[i][2];
    A->ha[i+1][3]=A->ha[i][3]+A->ha[i][5];
    A->ha[i][2]=A->ha[i][0];
    A->ha[i][5]=A->ha[i][3];
  } /* for i */
  A->ha[A->N-1][2]=A->ha[A->N-1][0];
  A->ha[A->N-1][5]=A->ha[A-N-1][3];
  i=A->rnr[A->NZ-1];
  j=A->cnr[A->NZ-1];
  xp=A->data[A->NZ-1];
  A->rnr[A->NZ-1]=-1;
  A->k=A->N;
  for(i3=1; i3<A->NZ; i3++) {
    i1=A->ha[i][2]+1;
    A->ha[i][2]=i1;
    i=A->rnr[i1];
    A->rnr[i1]=-1;
    z=A->data[i1];
    A->data[i1]=xp;
    xp=z;
    j1=A->cnr[i1];
    A->cnr[i1]=j;
    j=j1;
    if (i<=0) {
      do {
	i2=A->ha[k][0];
	A->k--;
	i=A->rnr[i2];
      } while (i<0);
      A->rnr[i2]=-1;
      xp=A->data[i2];
      j=A->cnr[i2];
    } /* if */
  } /* for i3 */
  i1=A->ha[i][2]+1;
  A->ha[i][2]=i1;
  A->data[i1]=xp;
  A->cnr[i1]=j;
  for(i=1; i<A->N; i++) {
    j1=A->ha[i][0]+1;
    j2=A->ha[i][2];
    for(j3=j1; j3<j2; j3++) {
      j=A->cnr[j3];
      k=A->ha[j][5]+1;
      A->rnr[k]=i;
      A->ha[j][5]=k;
    } /* for j3 */
  } /* for i */
} /* SPorder */


