#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpi.h>

#define min(x1,x2) ((x1) > (x2) ? (x2):(x1))

void read_xv8_(char*,int*,double*,double*,double*,int);
double findshort(double a, double b, double *lr, double BOXSIZE);

#define PI 3.141592653589
#define ALPHA 1.0/137.035999
#define HBARC 197.32693

int main(int argc, char *argv[]) {

  MPI_Status status;
  int PART, binnum, n, j, k, l, iter, flag, ionnum[5], dt, dumpi, ierr, rank, count, compnum;
  long long startt, stopt, confnum;
  double rho, BOXSIZE, dump, relpos[3], r, bind, time8;
  char line[300], filename[100];
  FILE *out[15], *ionfile;

  if(argc != 9) {
      fprintf(stderr,"Syntax Error: %s <ion.dat file> <output base> <# of particles> <density> <# of bins> <start time> <end time> <data cadence>\n",argv[0]);
      return 0;
    }

  ierr=MPI_Init(&argc,&argv);
  ierr=MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  ierr=MPI_Comm_size(MPI_COMM_WORLD,&count);

  ionfile=fopen(argv[1],"r");
  PART=atoi(argv[3]);
  sscanf(argv[4],"%lf",&rho);
  binnum=atoi(argv[5]);
  startt=atol(argv[6]);
  stopt=atol(argv[7]);
  dt=atoi(argv[8]);
  compnum=PART/count;

  confnum=(stopt-startt)/dt;

  BOXSIZE=pow(PART/rho,1.0/3.0);
  bind=BOXSIZE/(binnum*2.0);

  fgets(line,300,ionfile);

  double charge[PART], mass[PART];

  for(j=0;j<5;j++)
    ionnum[j]=0;

  for(j=0;j<PART;j++) {
    fgets(line,300,ionfile);
    sscanf(line,"%d %lf %lf",&dumpi,&charge[j],&mass[j]);
    if(charge[j]<3.0)
      ionnum[0]++;
    else if(charge[j]<7.0)
      ionnum[1]++;
    else if(charge[j]<9.0)
      ionnum[2]++;
    else if(charge[j]<11.0)
      ionnum[3]++;
    else 
      ionnum[4]++;
  }

  double vel[3*PART], pos[3*PART];

  int bins[15][binnum], recv[15][binnum];

  for(l=0;l<binnum;l++) 
    for(j=0;j<15;j++)
      bins[j][l]=0;
  
  for(iter=0;iter<confnum;iter++) {
    
    if(rank==0) {
      sprintf(filename,"md.ckpt.%011ld",startt+dt*iter);
      printf("%s\n",filename);
      k=strlen(filename);
      read_xv8_(filename,&PART,&time8,&pos[0],&vel[0],k);
    }
    
    ierr=MPI_Bcast(pos,3*PART,MPI_DOUBLE,0,MPI_COMM_WORLD);

    for(j=rank*compnum;j<(rank+1)*compnum;j++) { 
      for(k=0;k<j;k++) {
	for(n=0;n<3;n++) {
	  relpos[n]=findshort(pos[3*j+n],pos[3*k+n],&dump,BOXSIZE);
	}
	r=sqrt(relpos[0]*relpos[0]+relpos[1]*relpos[1]+relpos[2]*relpos[2]);
	flag=0;
	for(n=0;n<binnum;n++) {
	  if(flag==0) {
	    if(r<(n+1)*bind) {
	      if(charge[j]<3.0) {
		if(charge[k]<3.0)
		  bins[0][n]++;
		else if(charge[k]<7.0) 
		  bins[1][n]++;
		else if(charge[k]<9.0) 
		  bins[2][n]++;
		else if(charge[k]<11.0) 
		  bins[3][n]++;
		else 
		  bins[4][n]++;
	      }
	      else if(charge[j]<7.0) {
		if(charge[k]<3.0)
		  bins[1][n]++;
		else if(charge[k]<7.0) 
		  bins[5][n]++;
		else if(charge[k]<9.0) 
		  bins[6][n]++;
		else if(charge[k]<11.0) 
		  bins[7][n]++;
		else 
		  bins[8][n]++;
	      }
	      else if(charge[j]<9.0) {
		if(charge[k]<3.0)
		  bins[2][n]++;
		else if(charge[k]<7.0) 
		  bins[6][n]++;
		else if(charge[k]<9.0) 
		  bins[9][n]++;
		else if(charge[k]<11.0) 
		  bins[10][n]++;
		else 
		  bins[11][n]++;
	      }
	      else if(charge[j]<11.0) {
		if(charge[k]<3.0)
		  bins[3][n]++;
		else if(charge[k]<7.0) 
		  bins[7][n]++;
		else if(charge[k]<9.0) 
		  bins[10][n]++;
		else if(charge[k]<11.0) 
		  bins[12][n]++;
		else 
		  bins[13][n]++;
	      }
	      else {
		if(charge[k]<3.0)
		  bins[4][n]++;
		else if(charge[k]<7.0) 
		  bins[8][n]++;
		else if(charge[k]<9.0) 
		  bins[11][n]++;
		else if(charge[k]<11.0) 
		  bins[13][n]++;
		else 
		  bins[14][n]++;
	      }
	      flag=1;
	    }
	  }
	}
      }
    }
  }
  
  ierr=MPI_Allreduce(bins,recv,15*binnum,MPI_INT,MPI_SUM,MPI_COMM_WORLD);  

  /* if(rank==0) { */
  /*   if(ionnum[4]>0) { */
  /*     for(n=binnum-2;n<binnum;n++) { */
  /* 	r=(double)n/(double)binnum*BOXSIZE/2.0; */
  /* 	printf("%lf %e\n",r,(double)recv[9][n]*PART/(2*PI*r*r*rho*bind*ionnum[2]*ionnum[2]*confnum)); */
  /*     } */
  /*     for(n=binnum-2;n<binnum;n++) { */
  /* 	r=(double)n/(double)binnum*BOXSIZE/2.0; */
  /* 	printf("%lf %e\n",r,(double)recv[11][n]*PART/(4*PI*r*r*rho*bind*ionnum[2]*ionnum[4]*confnum)); */
  /*     } */
  /*     for(n=binnum-2;n<binnum;n++) { */
  /* 	r=(double)n/(double)binnum*BOXSIZE/2.0; */
  /* 	printf("%lf %e\n",r,(double)recv[14][n]*PART/(2*PI*r*r*rho*bind*ionnum[4]*ionnum[4]*confnum)); */
  /*     } */
  /*     printf("%d %d\n",ionnum[2],ionnum[4]); */
  /*   } */
  /*   else { */
  /*     for(n=binnum-2;n<binnum;n++) { */
  /* 	r=(double)n/(double)binnum*BOXSIZE/2.0; */
  /* 	printf("%lf %e\n",r,(double)recv[5][n]*PART/(2*PI*r*r*rho*bind*ionnum[1]*ionnum[1]*confnum)); */
  /*     } */
  /*     for(n=binnum-2;n<binnum;n++) { */
  /* 	r=(double)n/(double)binnum*BOXSIZE/2.0; */
  /* 	printf("%lf %e\n",r,(double)recv[6][n]*PART/(4*PI*r*r*rho*bind*ionnum[1]*ionnum[2]*confnum)); */
  /*     } */
  /*     for(n=binnum-2;n<binnum;n++) { */
  /* 	r=(double)n/(double)binnum*BOXSIZE/2.0; */
  /* 	printf("%lf %e\n",r,(double)recv[7][n]*PART/(4*PI*r*r*rho*bind*ionnum[1]*ionnum[3]*confnum)); */
  /*     } */
  /*     for(n=binnum-2;n<binnum;n++) { */
  /* 	r=(double)n/(double)binnum*BOXSIZE/2.0; */
  /* 	printf("%lf %e\n",r,(double)recv[9][n]*PART/(2*PI*r*r*rho*bind*ionnum[2]*ionnum[2]*confnum)); */
  /*     } */
  /*     for(n=binnum-2;n<binnum;n++) { */
  /* 	r=(double)n/(double)binnum*BOXSIZE/2.0; */
  /* 	printf("%lf %e\n",r,(double)recv[10][n]*PART/(4*PI*r*r*rho*bind*ionnum[2]*ionnum[3]*confnum)); */
  /*     } */
  /*     for(n=binnum-2;n<binnum;n++) { */
  /* 	r=(double)n/(double)binnum*BOXSIZE/2.0; */
  /* 	printf("%lf %e\n",r,(double)recv[12][n]*PART/(2*PI*r*r*rho*bind*ionnum[3]*ionnum[3]*confnum)); */
  /*     } */
  /*     printf("%d %d %d\n",ionnum[1],ionnum[2],ionnum[3]); */
  /*   } */
  /* } */
  if(rank==0) {
    if(ionnum[0]*ionnum[0]>0) {
      sprintf(filename,"%s.He-He.out",argv[2]);
      out[1]=fopen(filename,"w");
      for(n=0;n<binnum;n++) {
	r=(double)n/(double)binnum*BOXSIZE/2.0;
	fprintf(out[1],"%lf %e\n",r,(double)recv[0][n]*PART/(2*PI*r*r*rho*bind*ionnum[0]*ionnum[0]*confnum));
      }
    }
    if(ionnum[0]*ionnum[1]>0) {
      sprintf(filename,"%s.He-C.out",argv[2]);
      out[0]=fopen(filename,"w");
      for(n=0;n<binnum;n++) {
	r=(double)n/(double)binnum*BOXSIZE/2.0;
	fprintf(out[0],"%lf %e\n",r,(double)recv[1][n]*PART/(4*PI*r*r*rho*bind*ionnum[0]*ionnum[1]*confnum));
      }
    }
    if(ionnum[0]*ionnum[2]>0) {
      sprintf(filename,"%s.He-O.out",argv[2]);
      out[0]=fopen(filename,"w");
      for(n=0;n<binnum;n++) {
	r=(double)n/(double)binnum*BOXSIZE/2.0;
	fprintf(out[0],"%lf %e\n",r,(double)recv[2][n]*PART/(4*PI*r*r*rho*bind*ionnum[0]*ionnum[2]*confnum));
      }
    }
    if(ionnum[0]*ionnum[3]>0) {
      sprintf(filename,"%s.He-Ne.out",argv[2]);
      out[0]=fopen(filename,"w");
      for(n=0;n<binnum;n++) {
	r=(double)n/(double)binnum*BOXSIZE/2.0;
	fprintf(out[0],"%lf %e\n",r,(double)recv[3][n]*PART/(4*PI*r*r*rho*bind*ionnum[0]*ionnum[3]*confnum));
      }
    }
    if(ionnum[0]*ionnum[4]>0) {
      sprintf(filename,"%s.He-Se.out",argv[2]);
      out[0]=fopen(filename,"w");
      for(n=0;n<binnum;n++) {
	r=(double)n/(double)binnum*BOXSIZE/2.0;
	fprintf(out[0],"%lf %e\n",r,(double)recv[4][n]*PART/(4*PI*r*r*rho*bind*ionnum[0]*ionnum[4]*confnum));
      }
    }
    if(ionnum[1]*ionnum[1]>0) {
      sprintf(filename,"%s.C-C.out",argv[2]);
      out[0]=fopen(filename,"w");
      for(n=0;n<binnum;n++) {
	r=(double)n/(double)binnum*BOXSIZE/2.0;
	fprintf(out[0],"%lf %e\n",r,(double)recv[5][n]*PART/(2*PI*r*r*rho*bind*ionnum[1]*ionnum[1]*confnum));
      }
    }
    if(ionnum[1]*ionnum[2]>0) {
      sprintf(filename,"%s.C-O.out",argv[2]);
      out[0]=fopen(filename,"w");
      for(n=0;n<binnum;n++) {
	r=(double)n/(double)binnum*BOXSIZE/2.0;
	fprintf(out[0],"%lf %e\n",r,(double)recv[6][n]*PART/(4*PI*r*r*rho*bind*ionnum[1]*ionnum[2]*confnum));
      }
    }
    if(ionnum[1]*ionnum[3]>0) {
      sprintf(filename,"%s.C-Ne.out",argv[2]);
      out[0]=fopen(filename,"w");
      for(n=0;n<binnum;n++) {
	r=(double)n/(double)binnum*BOXSIZE/2.0;
	fprintf(out[0],"%lf %e\n",r,(double)recv[7][n]*PART/(4*PI*r*r*rho*bind*ionnum[1]*ionnum[3]*confnum));
      }
    }
    if(ionnum[1]*ionnum[4]>0) {
      sprintf(filename,"%s.C-Se.out",argv[2]);
      out[0]=fopen(filename,"w");
      for(n=0;n<binnum;n++) {
	r=(double)n/(double)binnum*BOXSIZE/2.0;
	fprintf(out[0],"%lf %e\n",r,(double)recv[8][n]*PART/(4*PI*r*r*rho*bind*ionnum[1]*ionnum[4]*confnum));
      }
    }
    if(ionnum[2]*ionnum[2]>0) {
      sprintf(filename,"%s.O-O.out",argv[2]);
      out[0]=fopen(filename,"w");
      for(n=0;n<binnum;n++) {
	r=(double)n/(double)binnum*BOXSIZE/2.0;
	fprintf(out[0],"%lf %e\n",r,(double)recv[9][n]*PART/(2*PI*r*r*rho*bind*ionnum[2]*ionnum[2]*confnum));
      }
    }
    if(ionnum[2]*ionnum[3]>0) {
      sprintf(filename,"%s.O-Ne.out",argv[2]);
      out[0]=fopen(filename,"w");
      for(n=0;n<binnum;n++) {
	r=(double)n/(double)binnum*BOXSIZE/2.0;
	fprintf(out[0],"%lf %e\n",r,(double)recv[10][n]*PART/(4*PI*r*r*rho*bind*ionnum[2]*ionnum[3]*confnum));
      }
    }
    if(ionnum[2]*ionnum[4]>0) {
      sprintf(filename,"%s.O-Se.out",argv[2]);
      out[0]=fopen(filename,"w");
      for(n=0;n<binnum;n++) {
	r=(double)n/(double)binnum*BOXSIZE/2.0;
	fprintf(out[0],"%lf %e\n",r,(double)recv[11][n]*PART/(4*PI*r*r*rho*bind*ionnum[2]*ionnum[4]*confnum));
      }
    }
    if(ionnum[3]*ionnum[3]>0) {
      sprintf(filename,"%s.Ne-Ne.out",argv[2]);
      out[0]=fopen(filename,"w");
      for(n=0;n<binnum;n++) {
	r=(double)n/(double)binnum*BOXSIZE/2.0;
	fprintf(out[0],"%lf %e\n",r,(double)recv[12][n]*PART/(2*PI*r*r*rho*bind*ionnum[3]*ionnum[3]*confnum));
      }
    }
    if(ionnum[3]*ionnum[4]>0) {
      sprintf(filename,"%s.Ne-Se.out",argv[2]);
      out[0]=fopen(filename,"w");
      for(n=0;n<binnum;n++) {
	r=(double)n/(double)binnum*BOXSIZE/2.0;
	fprintf(out[0],"%lf %e\n",r,(double)recv[13][n]*PART/(4*PI*r*r*rho*bind*ionnum[3]*ionnum[4]*confnum));
      }
    }
    if(ionnum[4]*ionnum[4]>0) {
      sprintf(filename,"%s.Se-Se.out",argv[2]);
      out[0]=fopen(filename,"w");
      for(n=0;n<binnum;n++) {
	r=(double)n/(double)binnum*BOXSIZE/2.0;
	fprintf(out[0],"%lf %e\n",r,(double)recv[14][n]*PART/(2*PI*r*r*rho*bind*ionnum[4]*ionnum[4]*confnum));
      }
    }
  }
}

double findshort(double a, double b, double *lr, double BOXSIZE) {
  double small, sign;
  
  // Say the closest neighbor is in the left box
  small=fabsf(a-(b-BOXSIZE));
  sign=1.0;
  
  // Check to see if it's in the middle box
  if(fabsf(a-b)<=fabsf(small)) {
    small=fabsf(a-b);
    if(small!=0.0)
      sign=(a-b)/small;
    else
      sign=1.0;
  }
  
  // Finally, see if it's in the right box
  if(fabsf(a-(b+BOXSIZE))<=fabsf(small)) {
    small=fabsf(a-(b+BOXSIZE));
    sign=-1.0;
  }
  
  // Give direction away from nearest neighbor
  *lr=sign;
  
  // Return the smallest value
  return small;
}
