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

  int PART, binnum, k, j, n, flag, coords, confnum, conf, ionnum[5], dumpi, nslice;
  float time4, dumpf;
  double time8;
  double maxd, bind, r, x, y, z, fastd, startt, stopt, cadence, dt, sliced;
  char filename[100], line[300];
  FILE *in, *finalin, *ionfile, *out[5];

  if(argc != 13) {
    fprintf(stderr,"Syntax Error: %s <ion.dat file> <output base> <# of particles> <# of bins> <# of slices> <max distance> <r, x, y, or z> <start time> <stop time> <Delta t> <data cadence> <Delta r threshold>\n",argv[0]);
    return 0;
  }

  ionfile=fopen(argv[1],"r");
  PART=atoi(argv[3]);
  binnum=atoi(argv[4]);
  nslice=atoi(argv[5]);
  sscanf(argv[6],"%lf",&maxd);
  coords=atoi(argv[7]);
  sscanf(argv[8],"%lf",&startt);
  sscanf(argv[9],"%lf",&stopt);
  sscanf(argv[10],"%lf",&dt);
  sscanf(argv[11],"%lf",&cadence);
  sscanf(argv[12],"%lf",&fastd);

  double pos[2][3*PART], skip[PART][3];
  double charge[PART], mass[PART];
  int bin[5][nslice][binnum], slice[PART], sliceion[5][nslice];

  bind=maxd/binnum;
  sliced=maxd/nslice;
  
  confnum=(int)((stopt-startt-dt)/cadence)+1;
  
  double diffusion[5][nslice], stdev[5];

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

  for(j=0;j<5;j++)
    for(n=0;n<nslice;n++) {
      diffusion[j][n]=0.0;
      sliceion[j][n]=0;
    }

  for(n=0;n<binnum;n++) 
    for(k=0;k<5;k++)
      for(j=0;j<nslice;j++)
	bin[k][j][n]=0;

  for(j=0;j<PART;j++) 
    slice[j]=0;

  sprintf(filename,"md.ckpt.%011.0lf",startt);
  k=strlen(filename);
  read_xv8_(filename,&PART,&time8,&pos[0][0],&skip[0][0],k);


  for(j=0;j<PART;j++) {
    z=pos[0][3*j+2];
    flag=0;
    for(n=0;n<nslice;n++) {
      if(flag==0) {
	if(z<(n+1)*sliced) {
	  slice[j]=n;
	  flag=1;
	}
      }
    }
    if(charge[j]<3.0)
      sliceion[0][slice[j]]++;
    else if(charge[j]<7.0)
      sliceion[1][slice[j]]++;
    else if(charge[j]<9.0)
      sliceion[2][slice[j]]++;
    else if(charge[j]<11.0)
      sliceion[3][slice[j]]++;
    else 
      sliceion[4][slice[j]]++;
  }

  for(conf=0;conf<confnum;conf++) {
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
      flag=0;
      if(r>fastd) {
	if(charge[j]<3.0){
	  diffusion[0][slice[j]]+=r*r;
	}
	else if(charge[j]<7.0){
	  diffusion[1][slice[j]]+=r*r;
	}
	else if(charge[j]<9.0){
	  diffusion[2][slice[j]]+=r*r;
	}
	else if(charge[j]<11.0){
	  diffusion[3][slice[j]]+=r*r;
	}
	else{
	  diffusion[4][slice[j]]+=r*r;
	}
      }
      for(n=0;n<binnum;n++) {
	if(flag==0) {
	  if(r<(n+1)*bind) {
	    if(charge[j]<3.0)
	      bin[0][slice[j]][n]++;
	    else if(charge[j]<7.0)
	      bin[1][slice[j]][n]++;
	    else if(charge[j]<9.0)
	      bin[2][slice[j]][n]++;
	    else if(charge[j]<11.0)
	      bin[3][slice[j]][n]++;
	    else
	      bin[4][slice[j]][n]++;
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
      printf("%lf",n*bind+bind/2.0);
      for(j=0;j<nslice;j++) {
	printf(" %e",(double)bin[0][j][n]/confnum);
      }
      printf("\n");
    }
    for(j=0;j<nslice;j++) {
      diffusion[0][j]=diffusion[0][j]/(6.0*dt*sliceion[0][j]*confnum);
      fprintf(out[0],"%lf %e\n",j*sliced+sliced/2.0,diffusion[0][j]);
    }
  }  
  if(ionnum[1]>0) {
    sprintf(filename,"%s.C.out",argv[2]);
    out[1]=fopen(filename,"w");
    for(n=0;n<binnum;n++) {
      printf("%lf",n*bind+bind/2.0);
      for(j=0;j<nslice;j++) {
	printf(" %e",(double)bin[1][j][n]/confnum);
      }
      printf("\n");
    }
    for(j=0;j<nslice;j++) {
      diffusion[1][j]=diffusion[1][j]/(6.0*dt*sliceion[1][j]*confnum);
      fprintf(out[1],"%lf %e\n",j*sliced+sliced/2.0,diffusion[1][j]);
    }
  }  
  if(ionnum[2]>0) {
    sprintf(filename,"%s.O.out",argv[2]);
    out[2]=fopen(filename,"w");
    for(n=0;n<binnum;n++) {
      printf("%lf",n*bind+bind/2.0);
      for(j=0;j<nslice;j++) {
	printf(" %e",(double)bin[2][j][n]/confnum);
      }
      printf("\n");
    }
    for(j=0;j<nslice;j++) {
      diffusion[2][j]=diffusion[2][j]/(6.0*dt*sliceion[2][j]*confnum);
      fprintf(out[2],"%lf %e\n",j*sliced+sliced/2.0,diffusion[2][j]);
    }
  }  
  if(ionnum[3]>0) {
    sprintf(filename,"%s.Ne.out",argv[2]);
    out[3]=fopen(filename,"w");
    for(n=0;n<binnum;n++) {
      printf("%lf",n*bind+bind/2.0);
      for(j=0;j<nslice;j++) {
	printf(" %e",(double)bin[3][j][n]/confnum);
      }
      printf("\n");
    }
    for(j=0;j<nslice;j++) {
      diffusion[3][j]=diffusion[3][j]/(6.0*dt*sliceion[3][j]*confnum);
      fprintf(out[3],"%lf %e\n",j*sliced+sliced/2.0,diffusion[3][j]);
    }
  }  
  if(ionnum[4]>0) {
    sprintf(filename,"%s.Fe.out",argv[2]);
    out[4]=fopen(filename,"w");
    for(n=0;n<binnum;n++) {
      printf("%lf",n*bind+bind/2.0);
      for(j=0;j<nslice;j++) {
	printf(" %e",(double)bin[4][j][n]/confnum);
      }
      printf("\n");
    }
    for(j=0;j<nslice;j++) {
      diffusion[4][j]=diffusion[4][j]/(6.0*dt*sliceion[4][j]*confnum);
      fprintf(out[4],"%lf %e\n",j*sliced+sliced/2.0,diffusion[4][j]);
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
