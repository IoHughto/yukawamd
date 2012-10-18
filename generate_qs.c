#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int sort_pass(int N, int *val1, int *val2, int *val3, int *val4);
int count_degens(int N, int nconf, int *val1);

int main(int argc, char *argv[]) {
  int j, k, l, N, nconf, count, N2max, dcount, dmax;
  int *val1, *val2, *val3, *val4;

  if(argc != 2) {
    fprintf(stderr,"Syntax Error: %s <n_i max>\n",argv[0]);
    return 1;
  }

  N=atoi(argv[1]);
  nconf=(2*N+1)*(2*N+1)*(2*N+1);

  N2max=3*N*N;

  val1=malloc(nconf*sizeof(int));
  val2=malloc(nconf*sizeof(int));
  val3=malloc(nconf*sizeof(int));
  val4=malloc(nconf*sizeof(int));

  count=0;

  for(j=0;j<2*N+1;j++) 
    for(k=0;k<2*N+1;k++)
      for(l=0;l<2*N+1;l++) {
	val2[count]=j-N;
	val3[count]=k-N;
	val4[count]=l-N;
	val1[count]=val2[count]*val2[count]+val3[count]*val3[count]+val4[count]*val4[count];
	count++;
      }

  if(count != nconf) {
    fprintf(stderr,"count != nconf: That's not good.\n");
    return 2;
  }

  while(sort_pass(nconf,&val1[0],&val2[0],&val3[0],&val4[0]));

  /*  for(j=0;j<nconf;j++) {
    printf("%d %d %d %d\n",val1[j],val2[j],val3[j],val4[j]);
    }*/

  dmax=0;

  for(j=0;j<=N2max;j++) {
    dcount=count_degens(j,nconf,&val1[0]);
    if(dcount>dmax)
      dmax=dcount;
  }

  printf("%d %d\n",N2max,dmax);

  for(j=0;j<=N2max;j++) {
    dcount=count_degens(j,nconf,&val1[0]);
    if(dcount>0) {
      printf("%d %d\n",j,dcount);
      for(k=0;k<nconf;k++)
	if(val1[k]==j)
	  printf("% 5d % 5d % 5d\n",val2[k],val3[k],val4[k]);
      printf("\n");
    }     
  }

}

int sort_pass(int N, int *val1, int *val2, int *val3, int *val4) {
  
  int swap, j, flag;

  flag=0;

  for(j=1;j<N;j++) {
    if(val1[j-1]>val1[j]) {
      swap=val1[j];
      val1[j]=val1[j-1];
      val1[j-1]=swap;
      swap=val2[j];
      val2[j]=val2[j-1];
      val2[j-1]=swap;
      swap=val3[j];
      val3[j]=val3[j-1];
      val3[j-1]=swap;
      swap=val4[j];
      val4[j]=val4[j-1];
      val4[j-1]=swap;
      flag=1;
    }
  }

  return flag;

}

int count_degens(int N, int nconf, int *val1) {
  int count, j;

  count=0;

  for(j=0;j<nconf;j++)
    if(val1[j]==N)
      count++;

  return count;
}
