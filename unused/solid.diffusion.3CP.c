#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define min(x1,x2) ((x1) > (x2) ? (x2):(x1))

void open_xv4a_(char*,int*,float*,float*,float*,int);
void read_xv4a_(char*,int*,float*,float*,float*,int);
void close_xv4a_(char*,int*,float*,float*,float*,int);

void read_xv8_(char*,int*,double*,double*,double*,int);

int main(int argc, char *argv[]) {

  int PART, binnum, cadence, k,j, n, flag, coords;
  float time4, dumpf;
  double time8;
  double maxd, bind, r, x, y, z;
  char filename[100], line[300];
  FILE *in, *finalin, *chargefile, *out[3];

  if(argc != 9) {
    fprintf(stderr,"Syntax Error: %s <input file 1> <input file 2> <charge file> <output file> <# of particles> <# of bins> <max distance> <r, x, y, or z>\n",argv[0]);
    return 0;
  }

  chargefile=fopen(argv[3],"r");
  sprintf(filename,"%s.C",argv[4]);
  out[0]=fopen(filename,"w");
  sprintf(filename,"%s.O",argv[4]);
  out[1]=fopen(filename,"w");
  sprintf(filename,"%s.Ne",argv[4]);
  out[2]=fopen(filename,"w");
  PART=atoi(argv[5]);
  binnum=atoi(argv[6]);
  sscanf(argv[7],"%lf",&maxd);
  coords=atoi(argv[8]);

  double pos[2][3*PART], skip[PART][3];
  double charge[PART];
  int bin[3][binnum];

  bind=maxd/binnum;

  for(j=0;j<PART;j++) {
    fgets(line,300,chargefile);
    sscanf(line,"%lf",&charge[j]);
  }

  fclose(chargefile);

  for(n=0;n<binnum;n++) {
    bin[0][n]=0;
    bin[1][n]=0;
    bin[2][n]=0;
  }

  sprintf(filename,"%s",argv[1]);
  k=strlen(filename);
  read_xv8_(filename,&PART,&time8,&pos[0][0],&skip[0][0],k);
  /* open_xv4a_(filename,&PART,&time4,&skip[0][0],&skip[0][0],k); */
  /* read_xv4a_(filename,&PART,&time4,&pos[0][0],&skip[0][0],k); */
  /* for(j=1;j<confnum-1;j++) { */
  /*   read_xv4a_(filename,&PART,&time4,&skip[0][0],&skip[0][0],k); */
  /*   if((j+1)%1000==0) */
  /*     printf("Skipped configuration %d\n",j+1); */
  /* } */
  /* read_xv4a_(filename,&PART,&time4,&pos[1][0],&skip[0][0],k); */
  /* close_xv4a_(filename,&PART,&time4,&skip[0][0],&skip[0][0],k); */
  sprintf(filename,"%s",argv[2]);
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
  
  for(n=0;n<binnum;n++) {
    fprintf(out[0],"%lf %d\n",n*bind+bind/2.0,bin[0][n]);
    fprintf(out[1],"%lf %d\n",n*bind+bind/2.0,bin[1][n]);
    fprintf(out[2],"%lf %d\n",n*bind+bind/2.0,bin[2][n]);
  }

}
