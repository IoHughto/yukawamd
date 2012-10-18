#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mpi.h"

#define min(x1,x2) ((x1) > (x2) ? (x2):(x1))

void open_xv4a_(char*,int*,float*,float*,float*,int);
void read_xv4a_(char*,int*,float*,float*,float*,int);
void close_xv4a_(char*,int*,float*,float*,float*,int);

void read_xv8_(char*,int*,double*,double*,double*,int);

int main(int argc, char *argv[]) {

  MPI_Status status;
  int PART, binnum, cadence, k,j, n, flag, coords, startt, stopt, dt, confnum, conf, ierr, count, rank, size;
  float time4, dumpf;
  double time8;
  double maxd, bind, r, x, y, z
, diffusion[3];
  char filename[100], line[300];
  FILE *in, *finalin, *chargefile, *out, *fastout;

  if(argc != 11) {
    fprintf(stderr,"Syntax Error: %s <charge file> <output file> <# of particles> <# of bins> <max distance> <r, x, y, or z> <start time> <stop time> <Delta t> <data cadence>\n",argv[0]);
    return 0;
  }

  ierr=MPI_Init(&argc, &argv);
  ierr=MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  ierr=MPI_Comm_size(MPI_COMM_WORLD,&count);

  chargefile=fopen(argv[1],"r");
  PART=atoi(argv[3]);
  binnum=atoi(argv[4]);
  sscanf(argv[5],"%lf",&maxd);
  coords=atoi(argv[6]);
  startt=atoi(argv[7]);
  stopt=atoi(argv[8]);
  dt=atoi(argv[9]);
  cadence=atoi(argv[10]);

  double pos[2][3*PART], skip[PART][3], fast[PART];
  double charge[PART];
  int bin[3][binnum], recv[3][binnum];

  bind=maxd/(binnum*2.0);

  confnum=(stopt-startt-dt)/cadence+1;

  for(j=0;j<PART;j++) {
    fgets(line,300,chargefile);
    sscanf(line,"%lf",&charge[j]);
  }

  fclose(chargefile);

  diffusion[0]=0.0; 
  diffusion[1]=0.0; 
  diffusion[2]=0.0;

  for(n=0;n<binnum;n++) {
    bin[0][n]=0;
    bin[1][n]=0;
    bin[2][n]=0;
  }

  for(conf=0;conf<confnum;conf++) {
    //    printf("%d %d\n",conf,confnum);
    sprintf(filename,"md.ckpt.%011d",startt+conf*cadence);
    k=strlen(filename);
    read_xv8_(filename,&PART,&time8,&pos[0][0],&skip[0][0],k);
    sprintf(filename,"md.ckpt.%011d",startt+conf*cadence+dt);
    k=strlen(filename);
    read_xv8_(filename,&PART,&time8,&pos[1][0],&skip[0][0],k);
    
    for(j=0;j<PART;j++) {
      x=min(min(fabsf(pos[1][3*j]-pos[0][3*j]),fabsf(pos[1][3*j]-pos[0][3*j]-maxd)),fabsf(pos[1][3*j]-pos[0][3*j]+maxd));
      y=min(min(fabsf(pos[1][3*j+1]-pos[0][3*j+1]),fabsf(pos[1][3*j+1]-pos[0][3*j+1]-maxd)),fabsf(pos[1][3*j+1]-pos[0][3*j+1]+maxd));
      z=min(min(fabsf(pos[1][3*j+2]-pos[0][3*j+2]),fabsf(pos[1][3*j+2]-pos[0][3*j+2]-maxd)),fabsf(pos[1][3*j+2]-pos[0][3*j+2]+maxd));
      r=sqrt(x*x+y*y+z*z);
      if(coords==1)
	r=x;
      else if(coords==2)
	r=y;
      else if(coords==3)
	r=z;
      
      if(r>61.5)
	fast[j]=1;
      else
	fast[j]=0;
    }
    sprintf(filename,"fast.%011d.xyz",startt+cadence*conf);
    fastout=fopen(filename,"w");
    fprintf(fastout,"%d\n[name]\n",PART);
    for(j=0;j<PART;j++) {
      if(fast[j]==1)
	fprintf(fastout,"C %e %e %e\n",pos[0][3*j],pos[0][3*j+1],pos[0][3*j+2]);
      else
	fprintf(fastout,"O %e %e %e\n",pos[0][3*j],pos[0][3*j+1],pos[0][3*j+2]);
    }
    fclose(fastout);
    size=PART/count;
    for(j=rank*size;j<(rank+1)*size;j++) {
      if((j+1)%1000==0 && rank==3)
	printf("%lf%\n",(double)(size*(count-1)*(j-size*(count-1)))/(double)(size*size*(count-1)+size*(size-1)/2)*100.0);
      for(k=0;k<j;k++) {
	x=min(min(fabsf(pos[0][3*j]-pos[0][3*k]),fabsf(pos[0][3*j]-pos[0][3*k]-maxd)),fabsf(pos[0][3*j]-pos[0][3*k]+maxd));
	y=min(min(fabsf(pos[0][3*j+1]-pos[0][3*k+1]),fabsf(pos[0][3*j+1]-pos[0][3*k+1]-maxd)),fabsf(pos[0][3*j+1]-pos[0][3*k+1]+maxd));
	z=min(min(fabsf(pos[0][3*j+2]-pos[0][3*k+2]),fabsf(pos[0][3*j+2]-pos[0][3*k+2]-maxd)),fabsf(pos[0][3*j+2]-pos[0][3*k+2]+maxd));
	r=sqrt(x*x+y*y+z*z);
	flag=0;
	if(r<maxd/2.0) {
	  for(n=0;n<binnum;n++) {
	    if(flag==0) {
	      if(r<(n+1)*bind) {
		if(fast[j]==1 && fast[k]==1)
		  bin[0][n]++;
		else if(fast[j]==0 && fast[k]==0)
		  bin[1][n]++;
		else
		  bin[2][n]++;
		flag=1;
	      }
	    }
	  }
	}
      }
    }
    printf("Proc %d done.\n",rank);
    ierr=MPI_Allreduce(bin,recv,3*binnum,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    if(rank==0) {
      sprintf(filename,"fast.%011d.CC",startt+cadence*conf);
      out=fopen(filename,"w");
      for(n=0;n<binnum;n++) {
	fprintf(out,"%e %d\n",(double)n/(double)binnum*maxd*sqrt(3.0)/4.0,recv[0][n]);
      }
      fclose(out);
      sprintf(filename,"fast.%011d.OO",startt+cadence*conf);
      out=fopen(filename,"w");
      for(n=0;n<binnum;n++) {
	fprintf(out,"%e %d\n",(double)n/(double)binnum*maxd*sqrt(3.0)/4.0,recv[1][n]);
      }
      fclose(out);
      sprintf(filename,"fast.%011d.CO",startt+cadence*conf);
      out=fopen(filename,"w");
      for(n=0;n<binnum;n++) {
	fprintf(out,"%e %d\n",(double)n/(double)binnum*maxd*sqrt(3.0)/4.0,recv[2][n]);
      }
      fclose(out);
    }
  } 
}
