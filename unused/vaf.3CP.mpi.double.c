#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mpi.h"

#define min(x1,x2) ((x1) > (x2) ? (x2):(x1))

void open_xv8a_(char*,int*,double*,double*,double*,int);
void read_xv8a_(char*,int*,double*,double*,double*,int);
void close_xv8a_(char*,int*,double*,double*,double*,int);

int main(int argc, char *argv[]) {

  MPI_Status status;
  int conf, j, n, confnum, cadence, PART, confmax, fmc, other, end, num[3], skipnum, k, offset, ierr, rank, count, compnum;
  char filename[100], line[300];
  double time8;
  float dumpf;
  FILE *in, *Cout, *Oout, *Neout, *chargefile;

  if(argc != 10) {
    fprintf(stderr,"Syntax Error: %s <input file> <charge file> <output file base> <# of configurations> <fm/c between ckpts> <max number of configuration steps> <# of particles> <# of configs to skip> <fm/c offset>\n",argv[0]);
    return 0;
  }

  ierr=MPI_Init(&argc, &argv);
  ierr=MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  ierr=MPI_Comm_size(MPI_COMM_WORLD,&count);

  chargefile=fopen(argv[2],"r");
  confnum=atoi(argv[4]);
  cadence=atoi(argv[5]);
  confmax=atoi(argv[6]);
  PART=atoi(argv[7]);
  skipnum=atoi(argv[8]);
  offset=atoi(argv[9]);

  sprintf(filename,"%s.C.out",argv[3]);
  Cout=fopen(filename,"w");
  sprintf(filename,"%s.O.out",argv[3]);
  Oout=fopen(filename,"w");
  sprintf(filename,"%s.Ne.out",argv[3]);
  Neout=fopen(filename,"w");

  compnum=PART/count;

  double vel[3*PART], skip[PART][3];
  double recvbuf[confmax][3*compnum];
  double vaf[3][confmax], charge[PART], pass, recv[3][confmax];
  int load[confnum];

  num[0]=0;
  num[1]=0;
  num[2]=0;

  for(j=0;j<PART;j++) {
    fgets(line,300,chargefile);
    sscanf(line,"%lf",&charge[j]);
    if(charge[j]<7.0)
      num[0]++;
    else if(charge[j]<9.0)
      num[1]++;
    else
      num[2]++;
  }

  fclose(chargefile);

  for(j=0;j<confmax;j++){
    vaf[0][j]=0.0;
    vaf[1][j]=0.0;
    vaf[2][j]=0.0;
  }

  for(j=0;j<confnum;j++)
    load[j]=0;

  if(rank == 0) {
    sprintf(filename,"%s",argv[1]);
    k=strlen(filename);
    open_xv8a_(filename,&PART,&time8,&skip[0][0],&skip[0][0],k);
  }

  if(rank == 0) {
    for(conf=0;conf<skipnum;conf++) {
      read_xv8a_(filename,&PART,&time8,&skip[0][0],&skip[0][0],k);
      printf("Skipping configuration %d\n",conf+1);
    }
  }

  for(conf=0;conf<confnum;conf++) {
    if(rank == 0) 
      printf("Configuration %d %e %e %e\n",conf+1,vaf[0][0],vaf[1][0],vaf[2][0]);
    if(load[conf]==0) {
      if(rank == 0) 
	read_xv8a_(filename,&PART,&time8,&skip[0][0],&vel[0],k);
      load[conf]=1;
      ierr=MPI_Scatter(vel,3*compnum,MPI_DOUBLE,recvbuf[conf%confmax],3*compnum,MPI_DOUBLE,0,MPI_COMM_WORLD);
    }
    end=min(conf+confmax,confnum);
    for(other=conf;other<end;other++) {
      if(load[other]==0) {
	if(rank == 0) 
	  read_xv8a_(filename,&PART,&time8,&skip[0][0],&vel[0],k);
	load[other]=1;
	ierr=MPI_Scatter(vel,3*compnum,MPI_DOUBLE,recvbuf[other%confmax],3*compnum,MPI_DOUBLE,0,MPI_COMM_WORLD);
      }
      for(j=0;j<compnum;j++) {
	pass=recvbuf[conf%confmax][3*j]*recvbuf[other%confmax][3*j]+recvbuf[conf%confmax][3*j+1]*recvbuf[other%confmax][3*j+1]+recvbuf[conf%confmax][3*j+2]*recvbuf[other%confmax][3*j+2];
	if(charge[compnum*rank+j]<7.0)
	  vaf[0][other-conf]+=pass;
	else if(charge[compnum*rank+j]<9.0)
	  vaf[1][other-conf]+=pass;
	else
	  vaf[2][other-conf]+=pass;
      } 
    }
  }
  close_xv8a_(filename,&PART,&time8,&skip[0][0],&skip[0][0],k);
  
  ierr=MPI_Allreduce(vaf,recv,3*confmax,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  if(rank == 0) {
    for(j=0;j<confmax;j++) {
      recv[0][j]/=(num[0]);
      fprintf(Cout,"%d %e\n",j*cadence,recv[0][j]);
      recv[1][j]/=(num[1]);
      fprintf(Oout,"%d %e\n",j*cadence,recv[1][j]);
      recv[2][j]/=(num[2]);
      fprintf(Neout,"%d %e\n",j*cadence,recv[2][j]);
    }
  }
  
  return 0;
  
}
