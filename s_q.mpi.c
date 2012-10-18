#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>
#include "mpi.h"

void open_xv8a_(char*,int*,float*,float*,float*,int);
void read_xv8a_(char*,int*,float*,float*,float*,int);
void close_xv8a_(char*,int*,float*,float*,float*,int);

int main(int argc, char *argv[]) {

  int j, k, l, n, PART, nconf, conf, N2max[2], dmax, qval, dcount, flag, ierr, rank, count, dconf, qcad;
  double tstart, tstop, cad, rho, q_max, time8, BOXSIZE, r;
  FILE *q_values;
  char filename[100],line[100];
  MPI_Status status;

  if(argc != 9) {
    fprintf(stderr,"Syntax error: %s <rho> <N> <q_max> <q value file> <tstart> <tstop> <data cadence> <n^2 cadence>\n",argv[0]);
    return 1;
  }

  ierr=MPI_Init(&argc,&argv);
  ierr=MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  ierr=MPI_Comm_size(MPI_COMM_WORLD,&count);

  sscanf(argv[1],"%lf",&rho);
  PART=atoi(argv[2]);
  sscanf(argv[3],"%lf",&q_max);
  sscanf(argv[5],"%lf",&tstart);
  sscanf(argv[6],"%lf",&tstop);
  sscanf(argv[7],"%lf",&cad);
  qcad=atoi(argv[8]);
  q_values=fopen(argv[4],"r");

  BOXSIZE=pow(PART/rho,1.0/3.0);
  N2max[0]=(int)(pow(q_max*BOXSIZE/(2.0*3.14159),2.0))+1;

  nconf=(int)((tstop-tstart)/cad);

  if(nconf%count != 0) {
    fprintf(stderr,"%d configurations cannont be broken up evenly between %d processors.\n",nconf,count);
    return 4;
  }

  dconf=nconf/count;

  flag=0;

  if(rank==0) {
    fgets(line,100,q_values);
    sscanf(line,"%d %d",&N2max[1],&dmax);
  }

  ierr=MPI_Bcast(&N2max[1],1,MPI_INT,0,MPI_COMM_WORLD);
  ierr=MPI_Bcast(&dmax,1,MPI_INT,0,MPI_COMM_WORLD);

  if(N2max[0]>N2max[1]) {
    if(rank==0)
      fprintf(stderr,"q_max is larger than the available number of q values in %s\n",argv[4]);
    return 3;
  }
  
  int qs[3][N2max[1]][dmax], qdegen[N2max[1]];
  double pos[PART][3], vel[3][PART], sq[5][N2max[1]], qd[3][N2max[1]][dmax];

  flag=0;

  for(j=0;j<N2max[0];j++)
    qdegen[j]=0;

  if(rank==0) 
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
  
  ierr=MPI_Bcast(qs,3*N2max[1]*dmax,MPI_INT,0,MPI_COMM_WORLD);
  ierr=MPI_Bcast(qdegen,N2max[1],MPI_INT,0,MPI_COMM_WORLD);
  
  for(j=0;j<N2max[0];j++) 
    for(n=0;n<5;n++)
      sq[n][j]=0.0;
  
  for(conf=0;conf<dconf;conf++) {
    for(k=0;k<N2max[0];k++) 
      for(n=0;n<2;n++) 
	sq[n][k]=0.0;
    sprintf(filename,"md.ckpt.%011.0lf",tstart+(dconf*rank+conf)*cad);
    k=strlen(filename);
    read_xv8_(filename,&PART,&time8,&pos[0],&vel[0],k);
    /*
    for(j=0;j<PART;j++)
      for(n=0;n<3;n++)
	pos[j][n]-=BOXSIZE/2.0;
    */
    for(j=0;j<PART;j++) 
      for(k=0;k<N2max[0];k++) {
	if(k%qcad==0) {
	  for(l=0;l<qdegen[k];l++) {
	    r=qs[0][k][l]*pos[j][0]+qs[1][k][l]*pos[j][1]+qs[2][k][l]*pos[j][2];
	    r*=2*3.14159/BOXSIZE;
	    sq[0][k]+=cos(r);
	    sq[1][k]+=sin(r);
	  }
	}
      }
    
    for(k=0;k<N2max[0];k++) {
      if(k%qcad==0) {
	sq[2][k]+=sq[0][k]*sq[0][k]+sq[1][k]*sq[1][k];
	sq[3][k]+=sq[0][k];
	sq[4][k]+=sq[1][k];
      }
    }
  }
  
  double recv[5][N2max[1]];

  ierr=MPI_Allreduce(sq,recv,5*N2max[1],MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD); 

  if(rank==0)
    for(k=0;k<N2max[0];k++) {
      if(k%qcad==0) {
	for(n=0;n<5;n++) 
	  if(qdegen[k]>0)
	    recv[n][k]/=(nconf*qdegen[k]*PART);
	if(qdegen[k]>0)
	  //	printf("% e % e % e % e % e\n",sqrt(k)*2.0*3.14159/BOXSIZE,sq[2][k]-sq[3][k]*sq[3][k]-sq[4][k]*sq[4][k],sq[2][k],sq[3][k],sq[4][k]);
	  printf("%e %e\n",sqrt(k)*2.0*3.14159/BOXSIZE,recv[2][k]);
	//printf("% e % e % e % e\n",sqrt(k)*2.0*3.14159/BOXSIZE,recv[0][k]*recv[0][k]+recv[1][k]*recv[1][k],recv[0][k],recv[1][k]);
      }
    }
  
  return 0;
}
