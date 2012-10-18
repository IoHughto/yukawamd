#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define min(x1,x2) ((x1) > (x2) ? (x2):(x1))

void open_xv4a_(char*,int*,float*,float*,float*,int);
void read_xv4a_(char*,int*,float*,float*,float*,int);
void close_xv4a_(char*,int*,float*,float*,float*,int);
double findshort(double a, double b, double *lr, double BOXSIZE);

void read_xv8_(char*,int*,double*,double*,double*,int);

int main(int argc, char *argv[]) {

  int PART, binnum, k, j, n, flag, coords, confnum, conf, ionnum[5], dumpi;
  float time4, dumpf;
  double time8;
  double maxd, bind, r, x, y, z, diffusion[2], fastd, startt, stopt, cadence, dt, dfac, lowZ, highZ, meanZ;
  char filename[100], line[300];
  FILE *in, *finalin, *ionfile, *out[2], *tout[2];

  if(argc != 14) {
    fprintf(stderr,"Syntax Error: %s <ion.dat file> <output base> <# of particles> <low Z> <high Z> <# of bins> <max distance> <r, x, y, z, or rho> <start time> <stop time> <Delta t> <data cadence> <Delta r threshold>\n",argv[0]);
    return 0;
  }

  ionfile=fopen(argv[1],"r");
  PART=atoi(argv[3]);
  sscanf(argv[4],"%lf",&lowZ);
  sscanf(argv[5],"%lf",&highZ);
  binnum=atoi(argv[6]);
  sscanf(argv[7],"%lf",&maxd);
  coords=atoi(argv[8]);
  sscanf(argv[9],"%lf",&startt);
  sscanf(argv[10],"%lf",&stopt);
  sscanf(argv[11],"%lf",&dt);
  sscanf(argv[12],"%lf",&cadence);
  sscanf(argv[13],"%lf",&fastd);

  meanZ=(lowZ+highZ)/2.0;

  switch(coords) {
  case 0:
    dfac=6.0;
    break;
  case 1:
    dfac=2.0;
    break;
  case 2:
    dfac=2.0;
    break;
  case 3:
    dfac=2.0;
    break;
  case 4:
    dfac=4.0;
    break;
  }

  double pos[2][3*PART], skip[PART][3];
  double charge[PART], mass[PART];
  int bin[2][binnum];

  bind=maxd/binnum;
  
  confnum=(int)((stopt-startt-dt)/cadence)+1;
  
  double diff[2][confnum], diffmean[2], stdev[2];

  for(j=0;j<2;j++)
    ionnum[j]=0;

  fgets(line,300,ionfile);

  for(j=0;j<PART;j++) {
    fgets(line,300,ionfile);
    sscanf(line,"%d %lf %lf",&dumpi,&charge[j],&mass[j]);
    if(charge[j]<meanZ)
      ionnum[0]++;
    else
      ionnum[1]++;
  }

  diffusion[0]=0.0; 
  diffusion[1]=0.0; 

  for(n=0;n<binnum;n++) {
    bin[0][n]=0;
    bin[1][n]=0;
  }

  for(conf=0;conf<confnum;conf++) {
    for(j=0;j<2;j++)
      diff[j][conf]=0.0;
    //    printf("%d %d\n",conf,confnum);
    sprintf(filename,"md.ckpt.%011.0lf",startt+conf*cadence);
    k=strlen(filename);
    read_xv8_(filename,&PART,&time8,&pos[0][0],&skip[0][0],k);
    sprintf(filename,"md.ckpt.%011.0lf",startt+conf*cadence+dt);
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
      else if(coords==4)
	r=sqrt(x*x+y*y);
      flag=0;
      if(r>fastd) {
	if(charge[j]<meanZ){
	  diffusion[0]+=r*r;
	  diff[0][conf]+=r*r;
	}
	else {
	  diffusion[1]+=r*r;
	  diff[1][conf]+=r*r;
	}
      }
      for(n=0;n<binnum;n++) {
	if(flag==0) {
	  if(r<(n+1)*bind) {
	    if(charge[j]<meanZ)
	      bin[0][n]++;
	    else 
	      bin[1][n]++;
	    flag=1;
	  }
	}
      }
    }
  } 
  
  if(ionnum[0]>0) {
    sprintf(filename,"%s.%d.out",argv[2],(int)lowZ);
    out[0]=fopen(filename,"w");
    sprintf(filename,"%s.%d.time.out",argv[2],(int)lowZ);
    tout[0]=fopen(filename,"w");
    for(n=0;n<binnum;n++) {
      fprintf(out[0],"%lf %e\n",n*bind+bind/2.0,(double)bin[0][n]/confnum);
    }
    diffmean[0]=diffusion[0]/(dfac*dt*ionnum[0]*confnum);
    stdev[0]=0.0;
    for(j=0;j<confnum;j++) 
      stdev[0]+=(diffmean[0]-diff[0][j]/(dfac*dt*ionnum[0]))*(diffmean[0]-diff[0][j]/(dfac*dt*ionnum[0]));
    stdev[0]=sqrt(stdev[0]/(confnum-1));
    for(conf=0;conf<confnum;conf++)
      fprintf(tout[0],"%11.0lf %16.14e\n",startt+conf*cadence,diff[0][conf]/(dfac*dt*ionnum[0]));
  }
  if(ionnum[1]>0) {
    sprintf(filename,"%s.%d.out",argv[2],(int)highZ);
    out[1]=fopen(filename,"w");
    sprintf(filename,"%s.%d.time.out",argv[2],(int)highZ);
    tout[1]=fopen(filename,"w");
    for(n=0;n<binnum;n++) {
      fprintf(out[1],"%lf %e\n",n*bind+bind/2.0,(double)bin[1][n]/confnum);
    }
    diffmean[1]=diffusion[1]/(dfac*dt*ionnum[1]*confnum);
    stdev[1]=0.0;
    for(j=0;j<confnum;j++) 
      stdev[1]+=(diffmean[1]-diff[1][j]/(dfac*dt*ionnum[1]))*(diffmean[1]-diff[1][j]/(dfac*dt*ionnum[1]));
    stdev[1]=sqrt(stdev[1]/(confnum-1));
    for(conf=0;conf<confnum;conf++)
      fprintf(tout[1],"%11.0lf %16.14e\n",startt+conf*cadence,diff[1][conf]/(dfac*dt*ionnum[1]));
  }
  if(ionnum[0]>0) {
    printf("D_low = %e +/- %e\n",diffmean[0],stdev[0]);
  }
  if(ionnum[1]>0) {
    printf("D_high = %e +/- %e\n",diffmean[1],stdev[1]);
  }
  return 0;
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
