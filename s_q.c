#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>

void open_xv8a_(char*,int*,float*,float*,float*,int);
void read_xv8a_(char*,int*,float*,float*,float*,int);
void close_xv8a_(char*,int*,float*,float*,float*,int);

int main(int argc, char *argv[]) {

  int j, k, l, n, PART, nconf, conf, N2max[2], dmax, qval, dcount, flag;
  double tstart, tstop, cad, rho, q_max, time8, BOXSIZE, r;
  FILE *q_values;
  char filename[100],line[100];

  if(argc != 8) {
    fprintf(stderr,"Syntax error: %s <rho> <N> <q_max> <q value file> <tstart> <tstop> <data cadence>\n",argv[0]);
    return 1;
  }

  sscanf(argv[1],"%lf",&rho);
  PART=atoi(argv[2]);
  sscanf(argv[3],"%lf",&q_max);
  sscanf(argv[5],"%lf",&tstart);
  sscanf(argv[6],"%lf",&tstop);
  sscanf(argv[7],"%lf",&cad);
  q_values=fopen(argv[4],"r");

  BOXSIZE=pow(PART/rho,1.0/3.0);
  N2max[0]=(int)(pow(q_max*BOXSIZE/(2.0*3.14159),2.0))+1;

  fgets(line,100,q_values);
  sscanf(line,"%d %d",&N2max[1],&dmax);

  if(N2max[0]>N2max[1]) {
    fprintf(stderr,"q_max is larger than the available number of q values in %s\n",argv[4]);
    return 3;
  }

  int qs[3][N2max[1]][dmax], qdegen[N2max[1]];
  double pos[PART][3], vel[3][PART], sq[5][N2max[1]], qd[3][N2max[1]][dmax];

  flag=0;

  for(j=0;j<N2max[0];j++)
    qdegen[j]=0;

  while(flag==0) {
    fgets(line,100,q_values);
    sscanf(line,"%d %d",&qval,&dcount);
    qdegen[qval]=dcount;
    for(j=0;j<dcount;j++) {
      fgets(line,100,q_values);
      sscanf(line,"%d %d %d\n",&qs[0][qval][j],&qs[1][qval][j],&qs[2][qval][j]);
      for(n=0;n<3;n++)
	qd[n][qval][j]=2*3.14159/BOXSIZE*qs[n][qval][j];
    }
    fgets(line,100,q_values);
    if(qval+1>N2max[0])
      flag=1;
  }
  /*
  for(j=0;j<N2max[0];j++) {
    printf("%d %d\n",j,qdegen[j]);
    if(qdegen[j]!=0) {
      for(k=0;k<qdegen[j];k++)
	printf("% 5d % 5d % 5d\n",qs[0][j][k],qs[1][j][k],qs[2][j][k]);
    }
    printf("\n");
  }
  */
  
  for(j=0;j<N2max[0];j++) 
    for(n=0;n<5;n++)
      sq[n][j]=0.0;

  nconf=(int)((tstop-tstart)/cad)+1;

  for(conf=0;conf<nconf;conf++) {
    for(k=0;k<N2max[0];k++) 
      for(n=0;n<2;n++) 
	sq[n][k]=0.0;
    sprintf(filename,"md.ckpt.%011.0lf",tstart+conf*cad);
    k=strlen(filename);
    read_xv8_(filename,&PART,&time8,&pos[0],&vel[0],k);
    for(j=0;j<PART;j++) 
      for(k=0;k<N2max[0];k++) 
	for(l=0;l<qdegen[k];l++) {
	  r=qd[0][k][l]*pos[j][0]+qd[1][k][l]*pos[j][1]+qd[2][k][l]*pos[j][2];
	  sq[0][k]+=cos(r)/qdegen[k];
	  sq[1][k]+=sin(r)/qdegen[k];	
	}
    
    for(k=0;k<N2max[0];k++) {
      sq[2][k]+=sq[0][k]*sq[0][k]+sq[1][k]*sq[1][k];
      sq[3][k]+=sq[0][k];
      sq[4][k]+=sq[1][k];
    }
  }
  
  for(k=0;k<N2max[0];k++) {
    for(n=2;n<5;n++) 
      sq[n][k]/=nconf;
    if(qdegen[k]!=0)
      //      printf("% e % e % e % e % e\n",sqrt(k)*2.0*3.14159/BOXSIZE,sq[2][k]-sq[3][k]*sq[3][k]-sq[4][k]*sq[4][k],sq[2][k],sq[3][k],sq[4][k]);
      //      printf("%e %e\n",k*2.0*3.14159/BOXSIZE,sq[2][k]);
      printf("% e % e % e % e\n",sqrt(k)*2.0*3.14159/BOXSIZE,sq[3][k]*sq[3][k]+sq[4][k]*sq[4][k],sq[3][k],sq[4][k]);
  }
  
  return 0;
}
