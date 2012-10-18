#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define min(x1,x2) ((x1) > (x2) ? (x2):(x1))

void read_xv8_(char*,int*,double*,double*,double*,int);
double findshort(double a, double b, double *lr, double BOXSIZE);

#define PI 3.141592653589
#define ALPHA 1.0/137.035999
#define HBARC 197.32693

int main(int argc, char *argv[]) {

  int PART, binnum, n, j, k, l, iter, flag, ionnum[5], startt, stopt, dt, confnum, dumpi;
  double rho, BOXSIZE, dump, relpos[3], r, bind, time8;
  char line[300], filename[100];
  FILE *out[5], *ionfile;

  if(argc != 9) {
      fprintf(stderr,"Syntax Error: %s <ion.dat file> <output base> <# of particles> <density> <# of bins> <start time> <end time> <data cadence>\n",argv[0]);
      return 0;
    }

  ionfile=fopen(argv[1],"r");
  PART=atoi(argv[3]);
  sscanf(argv[4],"%lf",&rho);
  binnum=atoi(argv[5]);
  startt=atoi(argv[6]);
  stopt=atoi(argv[7]);
  dt=atoi(argv[8]);

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

  double vel[PART][3], pos[PART][3];

  int bins[5][binnum];

  for(l=0;l<binnum;l++) 
    for(j=0;j<5;j++)
      bins[j][l]=0;
  
  for(iter=0;iter<confnum;iter++) {
    
    sprintf(filename,"md.ckpt.%011d",startt+dt*iter);
    k=strlen(filename);
    read_xv8_(filename,&PART,&time8,&pos[0][0],&vel[0][0],k);
    
    for(j=0;j<PART;j++) { 
      if((j+1)%1024==0)
	printf("%d\n",j+1);
      for(k=0;k<j;k++) {
	for(n=0;n<3;n++) {
	  relpos[n]=findshort(pos[j][n],pos[k][n],&dump,BOXSIZE);
	}
	r=sqrt(relpos[0]*relpos[0]+relpos[1]*relpos[1]+relpos[2]*relpos[2]);
	flag=0;
	for(n=0;n<binnum;n++) {
	  if(flag==0) {
	    if(r<(n+1)*bind) {
	      if(charge[j]<3.0 && charge[j]==charge[k])
		bins[0][n]++;
	      else if(charge[j]<7.0 && charge[j]==charge[k])
		bins[1][n]++;
	      else if(charge[j]<9.0 && charge[j]==charge[k])
		bins[2][n]++;
	      else if(charge[j]<11.0 && charge[j]==charge[k])
		bins[3][n]++;
	      else
		bins[4][n]++;
	      flag=1;
	    }
	  }
	}
      }
    }
  }
  
  if(ionnum[0]>0) {
    sprintf(filename,"%s.He.out",argv[2]);
    out[0]=fopen(filename,"w");
    for(n=0;n<binnum;n++) 
      fprintf(out[0],"%lf %lf\n",(double)n/(double)binnum*BOXSIZE/2.0,(double)bins[0][n]/(ionnum[0]*ionnum[0]*rho*rho/2.0*PI*pow(BOXSIZE,3.0)/(3.0*pow((double)binnum,3.0))*(pow(n+1,3.0)-pow(n,3.0)))*1.1547);
  }
  if(ionnum[1]>0) {
    sprintf(filename,"%s.C.out",argv[2]);
    out[1]=fopen(filename,"w");
    for(n=0;n<binnum;n++) 
      fprintf(out[1],"%lf %lf\n",(double)n/(double)binnum*BOXSIZE/2.0,(double)bins[1][n]/(ionnum[1]*ionnum[1]*rho*rho/2.0*PI*pow(BOXSIZE,3.0)/(3.0*pow((double)binnum,3.0))*(pow(n+1,3.0)-pow(n,3.0)))*1.1547);
  }
  if(ionnum[2]>0) {
    sprintf(filename,"%s.O.out",argv[2]);
    out[2]=fopen(filename,"w");
    for(n=0;n<binnum;n++) 
      fprintf(out[2],"%lf %lf\n",(double)n/(double)binnum*BOXSIZE/2.0,(double)bins[2][n]/(ionnum[2]*ionnum[2]*rho*rho/2.0*PI*pow(BOXSIZE,3.0)/(3.0*pow((double)binnum,3.0))*(pow(n+1,3.0)-pow(n,3.0)))*1.1547);
  }
  if(ionnum[3]>0) {
    sprintf(filename,"%s.Ne.out",argv[2]);
    out[3]=fopen(filename,"w");
    for(n=0;n<binnum;n++) 
      fprintf(out[3],"%lf %lf\n",(double)n/(double)binnum*BOXSIZE/2.0,(double)bins[3][n]/(ionnum[3]*ionnum[3]*rho*rho/2.0*PI*pow(BOXSIZE,3.0)/(3.0*pow((double)binnum,3.0))*(pow(n+1,3.0)-pow(n,3.0)))*1.1547);
  }
  if(ionnum[4]>0) {
    sprintf(filename,"%s.Fe.out",argv[2]);
    out[4]=fopen(filename,"w");
    for(n=0;n<binnum;n++) 
      fprintf(out[4],"%lf %lf\n",(double)n/(double)binnum*BOXSIZE/2.0,(double)bins[4][n]/(ionnum[4]*ionnum[4]*rho*rho/2.0*PI*pow(BOXSIZE,3.0)/(3.0*pow((double)binnum,3.0))*(pow(n+1,3.0)-pow(n,3.0)))*1.1547);
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
