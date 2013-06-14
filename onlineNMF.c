/*Author: Vihari Piratla
          Suraj B Malode
  Date: 14th June 2013
  Algorithm: [reference]
  Licence: MIT
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "levmar-2.6/levmar.h"

#ifndef LM_DBL_PREC
#error Demo program assumes that levmar has been compiled with double precision, see LM_DBL_PREC!
#endif

/*The product should not exceed 10**7 for avoiding seg fault*/
#define MAX_ROW 1000
#define MAX_COL 1000

int *x;                   /*temperorary vector read*/
int** F,G;                /*F and G to be computed*/

/*Returns random numbers between -1 and 1*/
double rand11()
{
  int r = rand();
  double s = ((double)r)/RAND_MAX;
  return (2*s-1);
}

/*Computes frobeniouys loss for given xi, F and gi
  Should be as efficient as possible, as being repetitively called*/
void frob_loss(double *p, double *err, int m, int n/*d*/, void *data){
  /*We can put F and x inside a structure and pass it, but that 
    instead making them global variables is more straighforward*/
  int i=0,j=0;
  /*similar to i of ref*/
  int timestep = *(int*)data;
  j = timestep;

  for(i=0;i<n;i++){
    err[i] = 0;
    
  }

}

int main(int argv,char** argc)
/*Reads input from std in line by line simulating the online behaviour
  Assumed that the matrix is in column major form i.e.
  each line consists of the column
 */
{
  int X[MAX_ROW][MAX_COLUMN];
  int d,n;

  /*Read d,n*/

  int r;                   /*The number of clusters (k means)*/
  
  if(argv>1)
    r = atoi(argc[1]);
  else
    r = 10;

  F = (int**)malloc(sizeof(int*)*d);
  int i;
  for (i=0;i<d;i++)
    F[i] = (int*)malloc(sizeof(int)*r);

  G = (int**)malloc(sizeof(int*)*r);
  for (i=0;i<r;i++)
    F[i] = (int*)malloc(sizeof(int)*n);

  int *g;
  g = (int*)malloc(r*sizeof(int));

  x = (int*)malloc(d*sizeof(int));

  int t;                  /*time counter: number of columns read so far*/
  
  int j;
  
  /*Initialize F*/
  for(i=0;i<d;i++)
    for(j=0;j<n;j++)
      /*F should be non negative*/
      F[i][j] = rand11()+1;

  while(true){
    for(j=0;j<n;j++){
    
      for(i=0;i<d;i++){
	if(scanf("%d",&x[i]) == EOF){
	  fprintf(stderr,"EOF reached before consuming the whole of input\n");
	  exit(0);
	}
      }
      
      /*TODO: This is inefficient as contradicting spatial cache coherence,
       so instead better take the other way and transpose while using it.*/
      for(i=0;i<d;i++){
	X[i][j] = x[i];
      }
      
      /*minimize frobenious loss with 
       constraints:
       1. sum over gi should equal 1
       2. gi should be non negative
      */
      
    }
  }
}
