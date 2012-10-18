#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define min(x1,x2) ((x1) > (x2) ? (x2):(x1))

double findshort(double a, double b, double *lr, double BOXSIZE);
void read_xv8_(char*,int*,double*,double*,double*,int);

int main(int argc, char *argv[]) {

  int PART, binnum, cadence, k, j, n, flag, coords, startt, stopt, dt, confnum, conf, ionnum[5], dumpi, slices, sum, i;
  float dumpf;
  double time8;
  double maxd, bind, r, x, y, z, rho, BOXSIZE, sliced;
  char filename[100], line[300];
  FILE *in, *finalin, *ionfile, *out[5];

  if(argc != 10) {
    fprintf(stderr,"Syntax Error: %s <ion.dat file> <output base> <# of particles> <density> <# of bins> <# of slices> <start time> <stop time> <data cadence>\n",argv[0]);
    return 0;
  }

  ionfile=fopen(argv[1],"r");
  PART=atoi(argv[3]);
  sscanf(argv[4],"%lf",&rho);
  binnum=atoi(argv[5]);
  slices=atoi(argv[6]);
  startt=atoi(argv[7]);
  stopt=atoi(argv[8]);
  cadence=atoi(argv[9]);

  double pos[3*PART], skip[PART][3];
  double charge[PART], mass[PART];
  int comps[5][slices];
  int bin[5][binnum];

  bind=1.0/binnum;
  sliced=1.0/slices;
  BOXSIZE=pow(PART/rho,1.0/3.0);
  confnum=(stopt-startt)/cadence+1;

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

  for(n=0;n<binnum;n++) {
    bin[0][n]=0;
    bin[1][n]=0;
    bin[2][n]=0;
    bin[3][n]=0;
    bin[4][n]=0;
  }

  for(conf=0;conf<confnum;conf++) {
    //    printf("%d %d\n",conf,confnum);
    for(n=0;n<slices;n++) {
      comps[0][n]=0;
      comps[1][n]=0;
      comps[2][n]=0;
      comps[3][n]=0;
      comps[4][n]=0;
    }
    sprintf(filename,"md.ckpt.%011d",startt+conf*cadence);
    k=strlen(filename);
    read_xv8_(filename,&PART,&time8,&pos[0],&skip[0][0],k);
    for(j=0;j<PART;j++) {
      flag=0;
      for(n=0;n<slices;n++) {
	if(pos[3*j+2]/BOXSIZE<(n+1)*sliced) {
	  if(flag==0) {
	    if(charge[j]<3.0)
	      comps[0][n]++;
	    else if(charge[j]<7.0)
	      comps[1][n]++;
	    else if(charge[j]<9.0)
	      comps[2][n]++;
	    else if(charge[j]<11.0)
	      comps[3][n]++;
	    else
	      comps[4][n]++;
	    flag=1;
	  }
	}
      }
    }
    for(n=0;n<slices;n++) {
      sum=comps[0][n]+comps[1][n]+comps[2][n]+comps[3][n]+comps[4][n];
      for(j=0;j<5;j++) {
	flag=0;
	for(i=0;i<binnum;i++) {
	  if((double)comps[j][n]/(double)sum<(i+1)*bind)
	    if(flag==0){
	      bin[j][i]++;
	      flag=1;
	    }
	}
      }
    }
  } 
  
  if(ionnum[0]>0) {
    sprintf(filename,"%s.He.out",argv[2]);
    out[0]=fopen(filename,"w");
    for(n=0;n<binnum;n++) {
      fprintf(out[0],"%lf %e\n",n*bind+bind/2.0,(double)bin[0][n]/confnum);
    }
  }

  if(ionnum[1]>0) {
    sprintf(filename,"%s.C.out",argv[2]);
    out[1]=fopen(filename,"w");
    for(n=0;n<binnum;n++) {
      fprintf(out[1],"%lf %e\n",n*bind+bind/2.0,(double)bin[1][n]/confnum);
    }
  }

  if(ionnum[2]>0) {
    sprintf(filename,"%s.O.out",argv[2]);
    out[2]=fopen(filename,"w");
    for(n=0;n<binnum;n++) {
      fprintf(out[2],"%lf %e\n",n*bind+bind/2.0,(double)bin[2][n]/confnum);
    }
  }

  if(ionnum[3]>0) {
    sprintf(filename,"%s.Ne.out",argv[2]);
    out[3]=fopen(filename,"w");
    for(n=0;n<binnum;n++) {
      fprintf(out[3],"%lf %e\n",n*bind+bind/2.0,(double)bin[3][n]/confnum);
    }
  }

  if(ionnum[4]>0) {
    sprintf(filename,"%s.Fe.out",argv[2]);
    out[4]=fopen(filename,"w");
    for(n=0;n<binnum;n++) {
      fprintf(out[4],"%lf %e\n",n*bind+bind/2.0,(double)bin[4][n]/confnum);
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
