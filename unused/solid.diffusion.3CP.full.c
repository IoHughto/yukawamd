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

  int PART, binnum, cadence, k,j, n, flag, coords, startt, stopt, dt, confnum, conf;
  float time4, dumpf;
  double time8;
  double maxd, bind, r, x, y, z
, diffusion[3];
  char filename[100], line[300];
  FILE *in, *finalin, *chargefile, *out[3];

  if(argc != 11) {
    fprintf(stderr,"Syntax Error: %s <charge file> <output file> <# of particles> <# of bins> <max distance> <r, x, y, or z> <start time> <stop time> <Delta t> <data cadence>\n",argv[0]);
    return 0;
  }

  chargefile=fopen(argv[1],"r");
  sprintf(filename,"%s.C",argv[2]);
  out[0]=fopen(filename,"w");
  sprintf(filename,"%s.O",argv[2]);
  out[1]=fopen(filename,"w");
  sprintf(filename,"%s.Ne",argv[2]);
  out[2]=fopen(filename,"w");
  PART=atoi(argv[3]);
  binnum=atoi(argv[4]);
  sscanf(argv[5],"%lf",&maxd);
  coords=atoi(argv[6]);
  startt=atoi(argv[7]);
  stopt=atoi(argv[8]);
  dt=atoi(argv[9]);
  cadence=atoi(argv[10]);

  double pos[2][3*PART], skip[PART][3];
  double charge[PART];
  int bin[3][binnum];

  bind=maxd/binnum;

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
      flag=0;
      if(r>16.0) {
	if(charge[j]<7.0)
	  diffusion[0]+=r*r;
	else if(charge[j]<9.0)
	  diffusion[1]+=r*r;
	else
	  diffusion[2]+=r*r;
      }
      for(n=0;n<binnum;n++) {
	if(flag==0) {
	  if(r<(n+1)*bind) {
	    if(charge[j]<7.0)
	      bin[0][n]++;
	    else if(charge[j]<9.0)
	      bin[1][n]++;
	    else
	      bin[2][n]++;
	    flag=1;
	  }
	}
      }
    }
  } 
  for(n=0;n<binnum;n++) {
    fprintf(out[0],"%lf %e\n",n*bind+bind/2.0,(double)bin[0][n]/confnum);
    fprintf(out[1],"%lf %e\n",n*bind+bind/2.0,(double)bin[1][n]/confnum);
    fprintf(out[2],"%lf %e\n",n*bind+bind/2.0,(double)bin[2][n]/confnum);
  }

  printf("D_C = %e\nD_O = %e\nD_Ne = %e\n",diffusion[0]/(6.0*dt*PART*confnum),diffusion[1]/(6.0*dt*PART*confnum),diffusion[2]/(6.0*dt*PART*confnum));
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
