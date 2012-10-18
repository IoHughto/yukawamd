#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char *argv[]) {

  int j, n, binnum, PART, flag;
  double BOXSIZE, rho;
  char line[300];
  FILE *in, *chargefile, *Cout, *Oout, *Neout;

  if(argc != 5) {
    fprintf(stderr,"Syntax Error: %s <positions> <charges> <# of bins> <# of particles>\n",argv[0]);
    return 0;
  }

  in=fopen(argv[1],"r");
  chargefile=fopen(argv[2],"r");
  Cout=fopen("mu.C.out","w");
  Oout=fopen("mu.O.out","w");
  Neout=fopen("mu.Ne.out","w");
  
  binnum=atoi(argv[3]);
  PART=atoi(argv[4]);

  rho=0.0000718;

  BOXSIZE=pow(PART/rho,1.0/3.0);

  int bins[3][binnum];
  double pos[PART][3], charge[PART];

  for(n=0;n<PART;n++) {
    fgets(line,300,chargefile);
    sscanf(line,"%lf",&charge[n]);
  }

  for(j=0;j<binnum;j++) {
    bins[0][j]=0;
    bins[1][j]=0;
    bins[2][j]=0;
  }

  double dumpd[3], time;
  float dumpf;

  fread(&dumpf,sizeof(float),1,in);
  fread(&time,sizeof(double),1,in);
  fread(&dumpf,sizeof(float),1,in);
  fread(&dumpf,sizeof(float),1,in);
  for(j=0;j<PART;j++) {
    fread(pos[j],sizeof(double),3,in);
    fread(dumpd,sizeof(double),3,in);
    flag=0;
    for(n=1;n<=binnum;n++)
      if(pos[j][2]>n*BOXSIZE/binnum) 
	flag=n;
    if(charge[j]<7.0)
      bins[0][flag]++;
    else if (charge[j]<9.0)
      bins[1][flag]++;
    else
      bins[2][flag]++;
  }
  fread(&dumpf,sizeof(float),1,in);

  for(n=0;n<binnum;n++) {
    fprintf(Cout,"%lf %d\n",BOXSIZE/binnum*(n+0.5),bins[0][n]);
    fprintf(Oout,"%lf %d\n",BOXSIZE/binnum*(n+0.5),bins[1][n]);
    fprintf(Neout,"%lf %d\n",BOXSIZE/binnum*(n+0.5),bins[2][n]);
  }
}
