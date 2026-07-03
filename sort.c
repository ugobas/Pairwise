//#include "nrutil.h"
#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
#define M 7       // size of subarrays sorted by straight insertion
#define NSTACK 50 // required auxiliary storage
#define NR_END 1
/**/

unsigned long *lvector(long nl, long nh);
void free_lvector(unsigned long *v, long nl, long nh);

void Sort(int *i_rank, float *X, unsigned long n){
  /* Sorts an array arr[1..n] into descending numerical order using the
     Quicksort algorithm. n is input; arr is replaced on output by its
     sorted  rearrangement. */

  float *arr=malloc((n+1)*sizeof(float));
  {// Copy array
    float *xx=X, *aa=arr+1;
    for(int i=0; i<n; i++){*aa=-*xx; aa++; xx++; i_rank[i]=i;}
  }

  unsigned long i, ir=n, j, k, l=1, i_r;
  int jstack=0;
  float a, temp;
  unsigned long *istack=lvector(1,NSTACK);
  for (;;) { // Insertion sort when subarray small enough.
    if (ir-l < M) {
      for (j=l+1;j<=ir;j++) {
	a=arr[j]; i_r=i_rank[j-1];
	for (i=j-1;i>=l;i--) {
	  if (arr[i] <= a) break;
	  arr[i+1]=arr[i]; i_rank[i]=i_rank[i-1];
	}
	arr[i+1]=a; i_rank[i]=i_r;
      }
      if (jstack == 0) break;
      ir=istack[jstack--]; //Pop stack and begin a new round of partitioning.
      l=istack[jstack--];
    } else {
      k=(l+ir) >> 1;
      /* Choose median of left, center and right elements as
	 partitioning  element a.Also rearrange so that
	 a[l]≤a[l+1]≤a[ir]. */
      SWAP(arr[k],arr[l+1])
	SWAP(i_rank[k-1],i_rank[l])
	if (arr[l] > arr[ir]) {
	  SWAP(arr[l],arr[ir])
	    SWAP(i_rank[l-1],i_rank[ir-1])
	    }
      if (arr[l+1] > arr[ir]) {
	SWAP(arr[l+1],arr[ir])
	  SWAP(i_rank[l],i_rank[ir-1])
	  }
      if (arr[l] > arr[l+1]) {
	SWAP(arr[l],arr[l+1])
	  SWAP(i_rank[l-1],i_rank[l])
	  }
      i=l+1; // Initialize pointers for partitioning.
      j=ir;
      a=arr[l+1]; i_r=i_rank[l]; //Partitioning  element.
      for (;;) { // Beginning  of  innermost  loop.
	do i++; while (arr[i] < a); //Scan up to find element > a.
	do j--; while (arr[j] > a); //Scan down to find element< a.
	if (j < i) break; // Pointers crossed. Partitioning  complete.
	SWAP(arr[i],arr[j]); // Exchange  elements.
	SWAP(i_rank[i-1],i_rank[j-1])
      } // End  of  innermost  loop.
      arr[l+1]=arr[j]; // Insert partitioning element.
      i_rank[l]=i_rank[j-1];
      arr[j]=a; i_rank[j-1]=i_r;
      jstack += 2;
      /* Push pointers to larger subarray on stack,
	 process smaller subarray immediately. */
      if (jstack > NSTACK){printf("NSTACK too small in sort."); exit(8);}
      if (ir-i+1 >= j-l) {
	istack[jstack]=ir;
	istack[jstack-1]=i;
	ir=j-1;
      } else {
	istack[jstack]=j-1;
	istack[jstack-1]=l;
	l=i;
      }
    }
  }
  free_lvector(istack,1,NSTACK);
  free(arr);
}

unsigned long *lvector(long nl, long nh)
     /* allocate an unsigned long vector with subscript range v[nl..nh] */
{
  unsigned long *v;
  
  v=(unsigned long *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(long)));
  if (!v){printf("allocation failure in lvector()"); exit(8);}
  return v-nl+NR_END;
}

void free_lvector(unsigned long *v, long nl, long nh)
     /* free an unsigned long vector allocated with lvector() */
{
  free((char *)(v+nl-NR_END));
}
