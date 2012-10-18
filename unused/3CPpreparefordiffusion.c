#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define min(x1,x2) ((x1) > (x2) ? (x2):(x1))

int main(int argc, char *argv[]) {

  int conf, j, n, confnum, cadence, PART, confmax, fmc, other, end, num[3], skipnum;
  char filename[100], line[300];
  double othervel[3], cm[3][3];
  float dumpf;
  FILE *in, *Cout, *Oout, *Neout, *chargefile;

  if(argc != 9) {
    fprintf(stderr,"Syntax Error: %s <input file> <charge file> <output file base> <# of configurations> <fm/c between ckpts> <max number of configuration steps> <# of particles> <# of configs to skip>\n",argv[0]);
    return 0;
  }

  in=fopen(argv[1],"r");
  chargefile=fopen(argv[2],"r");
  confnum=atoi(argv[4]);
  cadence=atoi(argv[5]);
  confmax=atoi(argv[6]);
  PART=atoi(argv[7]);
  skipnum=atoi(argv[8]);

  sprintf(filename,"%s.C.out",argv[3]);
  Cout=fopen(filename,"w");
  sprintf(filename,"%s.O.out",argv[3]);
  Oout=fopen(filename,"w");
  sprintf(filename,"%s.Ne.out",argv[3]);
  Neout=fopen(filename,"w");

  double vel[confmax][PART][3], vaf[3][confmax], charge[PART], skip[PART*3];
  int load[confnum];

  num[0]=0;
  num[1]=0;
  num[2]=0;

  for(j=0;j<3;j++) 
    for(n=0;n<3;n++)
      cm[j][n]=0.0;

  for(j=0;j<PART;j++) {
    fgets(line,300,chargefile);
    sscanf(line,"%lf",&charge[j]);
    if(charge[j]<7.0)
      num[0]++;
    else if(charge[j]<9.0)
      num[1]++;
    else
      num[2]++;
  }

  fclose(chargefile);

  for(j=0;j<confmax;j++){
    vaf[0][j]=0.0;
    vaf[1][j]=0.0;
    vaf[2][j]=0.0;
  }

  for(j=0;j<confnum;j++)
    load[j]=0;

  for(conf=0;conf<skipnum;conf++) {
    fread(skip,sizeof(double),3*PART,in);
    printf("Skipping configuration %d\n",conf+1);
  }

  for(conf=0;conf<confnum;conf++) {
    printf("Configuration %d\n",conf+1);
    if(load[conf]==0) {
      for(j=0;j<PART;j++) {
	fread(vel[conf%confmax][j],sizeof(double),3,in);
	for(n=0;n<3;n++) {
	if(charge[j]<7.0)
	  cm[0][n]+=vel[conf%confmax][j][n];
	else if(charge[j]<9.0)
	  cm[1][n]+=vel[conf%confmax][j][n];
	else
	  cm[2][n]+=vel[conf%confmax][j][n];
	}
	for(n=0;n<3;n++) {
	  cm[0][n]/=num[0];
	  cm[1][n]/=num[1];
	  cm[2][n]/=num[2];
	}
      }
      /* { */
      /* 	fgets(line,300,in); */
      /* 	sscanf(line,"%lf %lf %lf",&vel[conf%confmax][j][0],&vel[conf%confmax][j][1],&vel[conf%confmax][j][2]); */
      /* } */
      load[conf]=1;
    }
    end=min(conf+confmax,confnum);
    for(other=conf;other<end;other++) {
      if(load[other]==0) {
	for(j=0;j<PART;j++)

	  fread(vel[other%confmax][j],sizeof(double),3,in);
	/* { */
	/*   fgets(line,300,in); */
	/*   sscanf(line,"%lf %lf %lf",&vel[other%confmax][j][0],&vel[other%confmax][j][1],&vel[other%confmax][j][2]); */
	/* } */
	load[other]=1;
      }
      for(j=0;j<PART;j++) {
	if(charge[j]<7.0)
	  vaf[0][other-conf]+=(vel[conf%confmax][j][0]-cm[0][0])*(vel[other%confmax][j][0]-cm[0][0])+(vel[conf%confmax][j][1]-cm[0][1])*(vel[other%confmax][j][1]-cm[0][1])+(vel[conf%confmax][j][2]-cm[0][2])*(vel[other%confmax][j][2]-cm[0][2]);
	else if(charge[j]<9.0)
	  vaf[1][other-conf]+=(vel[conf%confmax][j][0]-cm[1][0])*(vel[other%confmax][j][0]-cm[1][0])+(vel[conf%confmax][j][1]-cm[1][1])*(vel[other%confmax][j][1]-cm[1][1])+(vel[conf%confmax][j][2]-cm[1][2])*(vel[other%confmax][j][2]-cm[1][2]);
	else 
	  vaf[2][other-conf]+=(vel[conf%confmax][j][0]-cm[2][0])*(vel[other%confmax][j][0]-cm[2][0])+(vel[conf%confmax][j][1]-cm[2][1])*(vel[other%confmax][j][1]-cm[2][1])+(vel[conf%confmax][j][2]-cm[2][2])*(vel[other%confmax][j][2]-cm[2][2]);
      }
    }
  }

  for(j=0;j<confmax;j++) {
    vaf[0][j]/=(num[0]);
    fprintf(Cout,"%d %e\n",j*cadence,vaf[0][j]);
    vaf[1][j]/=(num[1]);
    fprintf(Oout,"%d %e\n",j*cadence,vaf[1][j]);
    vaf[2][j]/=(num[2]);
    fprintf(Neout,"%d %e\n",j*cadence,vaf[2][j]);
  }

}
