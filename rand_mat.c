#include <stdio.h>
#include <stdlib.h>

double randnorm()
{
  int r = rand();
  double s = ((double)r)/RAND_MAX;
  return (s);
}

int main(int argc,char** argv){
  if(argc<2)
    exit(0);
  srand(time(NULL));
  int d = atoi(argv[1]);
  int n = atoi(argv[2]);
  
  printf("%d %d\n",d,n);
  int i,j;
  for(i=0;i<d;i++){
    for(j=0;j<n;j++)
      printf("%lf ",randnorm()*100);
    printf("\n");
  }
}
