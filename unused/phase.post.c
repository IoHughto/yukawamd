// Created by Joe Hughto
//
// The goal is to do all of the post-processing for phase diagram runs
//
// Routines:
// 1) g(r) calculation in slices
// 2) Diffusion for each direction, and r
// 3) Probability of slice with given composition vs. ion fraction

#include <stdio.h>
#include <math.h>
#include <string.h>

#define PI 3.141592653589
#define ALPHA 1.0/137.035999
#define HBARC 197.32693

void read_xv8_(char*,int*,double*,double*,double*,int);
double findshort(double a, double b, double *lr, double BOXSIZE);

int main(int argc, char *argv[]) {
  int routine[3], PART, binnum, slices, startt, stopt, cadence, confnum, ionnum[5], j, dumpi, k, iter, l, n, bflag, sflag;
  double dt, rho, bind, sliced, BOXSIZE, time8, rp[3], dumpd, r, pass;
  char line[300], filename[100];
  FILE *ionfile, *out[5];

  if(argc != 14) {
    fprintf(stderr,"Syntax Error: %s <ion.dat file> <output base> <# of particles> <density> <# of slices> <# of bins> <start time> <stop time> <data cadence> <g(r) calculation?> <diffusion?> <diffusion dt> <ion fraction probability?>\n",argv[0]);
    return 1;
  }

  ionfile=fopen(argv[1],"r");
  PART=atoi(argv[3]);
  sscanf(argv[4],"%lf",&rho); 
  slices=atoi(argv[5]); 
  binnum=atoi(argv[6]);
  startt=atoi(argv[7]);
  stopt=atoi(argv[8]);
  cadence=atoi(argv[9]);
  routine[0]=atoi(argv[10]);
  routine[1]=atoi(argv[11]);
  sscanf(argv[12],"%lf",&dt);
  routine[2]=atoi(argv[13]);

  double pos[PART][3], skip[PART][3], vel[PART][3];
  double charge[PART], mass[PART];
  int bin[5][slices][binnum];
  double binf[5][slices][binnum];

  BOXSIZE=pow(PART/rho,1.0/3.0);

  for(j=0;j<5;j++)
    ionnum[j]=0;

  fgets(line,300,ionfile);
    
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

  // g(r)
  if(routine[0]==1) {
    sliced=BOXSIZE/slices;
    bind=sliced/binnum;
    confnum=(stopt-startt)/cadence+1;
    for(j=0;j<5;j++) 
      for(k=0;k<slices;k++)
	for(l=0;l<binnum;l++)
	  bin[j][k][l]=0;
    for(iter=0;iter<confnum;iter++) {
      sprintf(filename,"md.ckpt.%011d",startt+cadence*iter);
      k=strlen(filename);
      read_xv8_(filename,&PART,&time8,&pos[0][0],&vel[0][0],k);

      for(j=0;j<PART;j++) {
	sflag=0;
	for(n=0;n<slices;n++) {
	  if(sflag==0) {
	    if(pos[j][2]<(n+1)*sliced) {
	      l=n;
	      sflag=1;
	    }
	  }
	}
	for(k=0;k<j;k++) {
	  for(n=0;n<3;n++)
	    rp[n]=findshort(pos[j][n],pos[k][n],&dumpd,BOXSIZE);
	  if(rp[2]<sliced) {
	    r=sqrt(rp[0]*rp[0]+rp[1]*rp[1]+rp[2]*rp[2]);
	    bflag=0;
	    for(n=0;n<binnum;n++){
	      if(bflag==0) {
		if(r<(n+1)*bind) {
		  if(charge[j]<3.0 && charge[j]==charge[k])
		    bin[0][l][n]++;
		  else if(charge[j]<7.0 && charge[j]==charge[k])
		    bin[1][l][n]++;
		  else if(charge[j]<9.0 && charge[j]==charge[k])
		    bin[2][l][n]++;
		  else if(charge[j]<11.0 && charge[j]==charge[k])
		    bin[3][l][n]++;
		  else
		    bin[4][l][n]++;
		  bflag=1;
		}
	      }
	    }
	  }
	}
      }
    }
    if(ionnum[0]>0) {
      sprintf(filename,"%s.gofr.He.out",argv[2]);
      out[0]=fopen(filename,"w");
      for(n=0;n<binnum;n++) {
	fprintf(out[0],"%lf ",(double)n/bind);
	for(l=0;l<slices;l++)
	  fprintf(out[0],"%lf ",(double)bin[0][l][n]/(ionnum[0]*ionnum[0]*rho*rho/2.0*PI*pow(BOXSIZE,3.0)/(3.0*pow((double)binnum,3.0))*(pow(n+1,3.0)-pow(n,3.0)))*1.1547);
	fprintf(out[0],"\n");
      }
    } 
    if(ionnum[1]>0) {
      sprintf(filename,"%s.gofr.C.out",argv[2]);
      out[1]=fopen(filename,"w");
      for(n=0;n<binnum;n++) {
	fprintf(out[1],"%e ",(double)n/bind);
	for(l=0;l<slices;l++)
	  fprintf(out[1],"%e ",(double)bin[1][l][n]/(ionnum[1]*ionnum[1]*rho*rho/2.0*PI*pow(BOXSIZE,3.0)/(3.0*pow((double)binnum,3.0))*(pow(n+1,3.0)-pow(n,3.0)))*1.1547);
	fprintf(out[1],"\n");
      }
    }
    if(ionnum[2]>0) {
      sprintf(filename,"%s.gofr.O.out",argv[2]);
      out[2]=fopen(filename,"w");
      for(n=0;n<binnum;n++) {
	fprintf(out[2],"%e ",(double)n/bind);
	for(l=0;l<slices;l++)
	  fprintf(out[2],"%e ",(double)bin[2][l][n]/(ionnum[2]*ionnum[2]*rho*rho/2.0*PI*pow(BOXSIZE,3.0)/(3.0*pow((double)binnum,3.0))*(pow(n+1,3.0)-pow(n,3.0)))*1.1547);
	fprintf(out[2],"\n");
      }
    }
    if(ionnum[3]>0) {
      sprintf(filename,"%s.gofr.Ne.out",argv[2]);
      out[3]=fopen(filename,"w");
      for(n=0;n<binnum;n++) {
	fprintf(out[3],"%e ",(double)n/bind);
	for(l=0;l<slices;l++)
	  fprintf(out[3],"%e ",(double)bin[3][l][n]/(ionnum[3]*ionnum[3]*rho*rho/2.0*PI*pow(BOXSIZE,3.0)/(3.0*pow((double)binnum,3.0))*(pow(n+1,3.0)-pow(n,3.0)))*1.1547);
	fprintf(out[3],"\n");
      }
    }
    if(ionnum[4]>0) {
      sprintf(filename,"%s.gofr.Fe.out",argv[2]);
      out[4]=fopen(filename,"w");
      for(n=0;n<binnum;n++) {
	fprintf(out[4],"%e ",(double)n/bind);
	for(l=0;l<slices;l++)
	  fprintf(out[4],"%e ",(double)bin[4][l][n]/(ionnum[4]*ionnum[4]*rho*rho/2.0*PI*pow(BOXSIZE,3.0)/(3.0*pow((double)binnum,3.0))*(pow(n+1,3.0)-pow(n,3.0)))*1.1547);
	fprintf(out[4],"\n");
      }
    }
  }
    
  // Diffusion
  if(routine[1]==1) {

  }

  // Slice composition probability
  if(routine[2]==1) {

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
