#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void read_xv8_(char*,int*,double*,double*,double*,int);

int main(int argc, char *argv[]) {
  int PART, binnum, k, j, transi, dumpi, confnum, conf, rect;
  double time8, trans, BOXSIZE, rho, startt, stopt, cadence;
  char line[300], filename[100];
  FILE *ionfile;

  if(argc != 8) {
    fprintf(stderr,"Syntax Error: %s <ion file> <density in 1/fm^3> <# of slices> <start time> <stop time> <data cadence> <rectangular?>\n",argv[0]);
    return 0;
  }
  
  ionfile=fopen(argv[1],"r");
  sscanf(argv[2],"%lf",&rho);
  binnum=atoi(argv[3]);
  sscanf(argv[4],"%lf",&startt);
  sscanf(argv[5],"%lf",&stopt);
  sscanf(argv[6],"%lf",&cadence);
  rect=atoi(argv[7]);

  fgets(line,300,ionfile);
  sscanf(line,"%d",&PART);

  int bins[3][binnum];
  double charge[PART], mass[PART], pos[PART][3], vel[PART][3];

  for(j=0;j<PART;j++) {
    fgets(line,300,ionfile);
    sscanf(line,"%d %lf %lf",&dumpi,&charge[j],&mass[j]);
  }

  confnum=(int)((stopt-startt)/cadence)+1;

  for(conf=0;conf<confnum;conf++) {

    for(j=0;j<binnum;j++) {
      bins[0][j]=0;
      bins[1][j]=0;
      bins[2][j]=0;
    }
    
    sprintf(filename,"md.ckpt.%011.0lf",startt+conf*cadence);
    k=strlen(filename);
    read_xv8_(filename,&PART,&time8,&pos[0][0],&vel[0][0],k);
    
    if(rect==0)
      BOXSIZE=pow(PART/rho,1.0/3.0);
    else
      BOXSIZE=2.0*pow(PART/(2.0*rho),1.0/3.0);

    for(j=0;j<PART;j++) {
      trans=pos[j][2]/BOXSIZE*binnum;
      transi=(int)floor(trans);
      if(charge[j]<7.0)
	bins[0][transi]++;
      else if(charge[j]<9.0)
	bins[1][transi]++;
      else
	bins[2][transi]++;
    }
    
    for(j=0;j<binnum;j++)
      printf("%d %d %d %d\n",bins[0][j],bins[1][j],bins[2][j],bins[0][j]+bins[1][j]+bins[2][j]);
    
  }
}
