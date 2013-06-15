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
#include <time.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "levmar-2.6/levmar.h"

#ifndef LM_DBL_PREC
#error This program assumes that levmar has been compiled with double precision, see LM_DBL_PREC!
#endif

/*The product should not exceed 10**7 for avoiding seg fault*/
#define MAX_ROW 1000
#define MAX_COL 1000
#define max(a,b) a>b?a:b
#define true 1

double *x;
/*F and G to be computed*/
gsl_matrix *F_gsl;
int n,d;                  /*n and d similar to ref*/
double opts[LM_OPTS_SZ], info[LM_INFO_SZ];

void write(FILE* f,gsl_matrix* m){
  int i1,i2;
  int b1=m->size1;
  int b2 = m->size2;
  for(i1=0;i1<b1;i1++){
    for(i2=0;i2<b2;i2++)
      fprintf(f,"%lf ",gsl_matrix_get(m,i1,i2));
    fprintf(f,"\n");
  }
}

/*writes in to R interpretable object*/
void write_R(FILE* f,gsl_matrix* m,const char* name){
  int i1,i2;
  int b1=m->size1;
  int b2 = m->size2;
  fprintf(f,"%s<-matrix(c(",name);
  for(i1=0;i1<b1;i1++){
    for(i2=0;i2<b2;i2++){
      if((i1!=(b1-1))||(i2!=(b2-1)))
	fprintf(f,"%lf ,",gsl_matrix_get(m,i1,i2));
      else
	fprintf(f,"%lf ",gsl_matrix_get(m,i1,i2));
    }
  }
  fprintf(f,"),nrow=%d,ncol=%d)\n",b1,b2);
}

/*Returns random numbers between -1 and 1*/
double rand11()
{
  int r = rand();
  double s = ((double)r)/RAND_MAX;
  return (2*s-1);
}

/*Computes frobeniouys loss for given xi, F and gi
  Should be as efficient as possible, as being repetitively called*/
void frob_loss(double *p, double *err, int m, int ds, void *data){
  register int i=0,j=0;
  gsl_matrix * F_gsl = (gsl_matrix*)data;
  fprintf(stderr,"\n");
  for(i=0;i<ds;i++){
    double temp = 0;
    for(j=0;j<m;j++){
      temp += gsl_matrix_get(F_gsl,i,j)*p[j];
    }
    err[i] = temp;
    fprintf(stderr,"%lf\n",temp);
  }
}

int main(int argv,char** argc)
/*Reads input from std in line by line simulating the online behaviour
  Assumed that the matrix is in column major form i.e.
  each line consists of the column
 */
{
  srand(time(NULL));
  int d,n;

  /*Read d,n*/
  scanf("%d %d",&d,&n);
  /*check for values*/
  if(d<n){
    fprintf(stderr,"d is supposed to be greater than n: \nfor the computation of g to be possible at any step\n");
    exit(0);
  }

  int r;                   /*The number of clusters (k means)*/
  
  if(argv>1)
    r = atoi(argc[1]);
  else{
    r = n;
    fprintf(stderr,"r value assumed to be %d\n",r);
  }

  gsl_matrix *G_gsl = gsl_matrix_alloc(r, n);
  gsl_matrix *X_gsl = gsl_matrix_alloc(d, n);
  int i;

  double *g;
  g = (double*)malloc(r*sizeof(double));
  x = (double*)malloc(d*sizeof(double));
  
  int t;                  /*time counter: number of columns read so far*/
  int j;

  gsl_matrix_view x_gsl = gsl_matrix_view_array(x, d, 1);
  gsl_matrix_view g_gsl = gsl_matrix_view_array(g, n, 1);
  /*velocity and hessian matrices*/  
  /*velocity is of dimension d*n and H is 
    of dimension n*n*/
  gsl_matrix *v_gsl = gsl_matrix_alloc(d, r);
  gsl_matrix *H_gsl = gsl_matrix_alloc(r, r);
  gsl_matrix *F_gsl = gsl_matrix_alloc(d, r);
  
  /*Initialize F*/
  for(i=0;i<d;i++)
    for(j=0;j<r;j++){
      /*F should be non negative*/
      double x = rand11()+1;
      gsl_matrix_set(F_gsl,i,j,x);
    }
  //printf("F:\n");
  //write(stdout,F_gsl);
  //printf("d: %d n: %d\n",d,n);

  while(true){
    for(j=0;j<n;j++){
    
      /*Under the assumptions that the data is in column major form*/
      for(i=0;i<d;i++){
	if(scanf("%lf",&x[i]) == EOF){
	  if(i!=0)
	    fprintf(stderr,"EOF reached before consuming the whole of input\n");
	  exit(0);
	}
      }
      
      for(i=0;i<d;i++)
	gsl_matrix_set(X_gsl,i,j, x[i]);
      
      /*minimize frobenious loss with 
       constraints:
       1. sum over gi should equal 1 [[for the moment this is neglected]]
       2. gi should be non negative
      */
      double lb[5], ub[5];
      
      for(i=0;i<n;i++){
	/*thereotically should be 0*/
	lb[i]=0.001;
	/*kind of infinity*/
	ub[i]=100000000.0;
      }
       
      for(i=0;i<r;i++){
	g[i] = rand11()+1;
      }
      //write(stderr,F_gsl);

      fprintf(stderr,"target:\n");
      for(i=0;i<d;i++)
	fprintf(stderr,"%lf ",x[i]);
      fprintf(stderr,"\n");
      
      //printf("main: n: %d d:%d\n",n,d);
      //int ret=dlevmar_dif(frob_loss, g, x, n, d, 1000, opts, info, NULL, NULL, (void*)F_gsl);
      int ret=dlevmar_bc_dif(frob_loss, g, x, n, d, lb, ub, NULL, 1000, opts, info, NULL, NULL, (void*)F_gsl);
      fprintf(stderr,"Converged in %d steps\n",ret);
      
      /*Updating of F*/
      /*Need some temp arrays for storing temporary results*/
           
      gsl_matrix* prod_F_H_gsl = gsl_matrix_alloc(d, r);
      gsl_matrix* grad_gsl = gsl_matrix_alloc(d, r);

      x_gsl = gsl_matrix_view_array(x, d, 1);
      g_gsl = gsl_matrix_view_array(g, n, 1);
  
      if(t==0){
	gsl_blas_dgemm (CblasNoTrans, CblasTrans,
			1.0, &x_gsl.matrix, &g_gsl.matrix,
			0.0, v_gsl);
	gsl_blas_dgemm (CblasNoTrans, CblasTrans,
			1.0, &g_gsl.matrix, &g_gsl.matrix,
			0.0, H_gsl);
      }
      else{
	/*The only difference with above module is the beta value*/
	gsl_blas_dgemm (CblasNoTrans, CblasTrans,
			1.0, &x_gsl.matrix, &g_gsl.matrix,
			1.0, v_gsl);
	gsl_blas_dgemm (CblasNoTrans, CblasTrans,
			1.0, &g_gsl.matrix, &g_gsl.matrix,
			1.0, H_gsl);
      }

      gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
		      1.0, F_gsl, H_gsl,
		      0.0, prod_F_H_gsl);

      /*grad = v-prod_F_H*/
      /*result stored in prod_F_H and bears negative result*/
      gsl_matrix_sub(prod_F_H_gsl,v_gsl);
	  
      gsl_blas_dtrsm(CblasRight,CblasLower,CblasNoTrans,CblasNonUnit,-1.0,H_gsl,prod_F_H_gsl);

      gsl_matrix_add(F_gsl,prod_F_H_gsl);
      int i1,i2;
      for(i1=0;i1<d;i1++)
	for(i2=0;i2<r;i2++)
	  gsl_matrix_set(F_gsl,i1,i2,max(0,gsl_matrix_get(F_gsl,i1,i2)));
      
      for(i1=0;i1<r;i1++){
	gsl_matrix_set(G_gsl,i1,j,g[i1]);
      }
    }
    /*for classic matrix dumping use write instead*/
    write_R(stdout,F_gsl,"F");
    write_R(stdout,G_gsl,"G");
    write_R(stdout,X_gsl,"X");
  }
}
