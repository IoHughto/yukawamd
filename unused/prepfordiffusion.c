#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define min(x1,x2) ((x1) > (x2) ? (x2):(x1))

int main(int argc, char *argv[]) {

  int conf, j, confnum, cadence, PART, confmax, fmc, other, end;
  char filename[100], line[300];
  double othervel[3];
  float dumpf;
  FILE *in, *out;

  if(argc != 7) {
    fprintf(stderr,"Syntax Error: %s <filename base> <output file> <# of configurations> <fm/c between ckpts> <max number of configuration steps> <# of particles>\n",argv[0]);
    return 0;
  }

  out=fopen(argv[2],"w");
  confnum=atoi(argv[3]);
  cadence=atoi(argv[4]);
  confmax=atoi(argv[5]);
  PART=atoi(argv[6]);

  double vel[confmax][PART][3], vaf[confmax];
  int load[confnum];

  for(j=0;j<confmax;j++)
    vaf[j]=0.0;

  for(j=0;j<confnum;j++)
    load[j]=0;

  for(conf=0;conf<confnum;conf++) {
    if(load[conf]==0) {
      fmc=(conf+1)*cadence;
      sprintf(filename,"%s.%011d",argv[1],fmc);
      printf("Loading %s\n",filename);
      in=fopen(filename,"r");
      for(j=0;j<PART;j++) {
	fgets(line,300,in);
	sscanf(line,"%lf %lf %lf",&vel[conf%confmax][j][0],&vel[conf%confmax][j][1],&vel[conf%confmax][j][2]);
      }
      fclose(in);
      load[conf]=1;
    }
    end=min(conf+confmax,confnum);
    for(other=conf;other<end;other++) {
      if(load[other]==0) {
	fmc=(other+1)*cadence;
	sprintf(filename,"%s.%011d",argv[1],fmc);
	printf("Loading %s\n",filename);
	in=fopen(filename,"r");
	for(j=0;j<PART;j++) {
	  fgets(line,300,in);
	  sscanf(line,"%lf %lf %lf",&vel[other%confmax][j][0],&vel[other%confmax][j][1],&vel[other%confmax][j][2]);
	}
	fclose(in);
	load[other]=1;
      }
      for(j=0;j<PART;j++) 
	vaf[other-conf]+=vel[conf%confmax][j][0]*vel[other%confmax][j][0]+vel[conf%confmax][j][1]*vel[other%confmax][j][1]+vel[conf%confmax][j][2]*vel[other%confmax][j][2];
    }
  }

  for(j=0;j<confmax;j++) {
    vaf[j]/=(PART*(confnum-j));
    fprintf(out,"%d %lf\n",j*cadence,vaf[j]);
  }

}
