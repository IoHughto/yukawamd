#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define min(x1,x2) ((x1) > (x2) ? (x2):(x1))

void open_xv4a_(char*,int*,double*,double*,double*,int);
void read_xv4a_(char*,int*,double*,double*,double*,int);
void close_xv4a_(char*,int*,double*,double*,double*,int);

int main(int argc, char *argv[]) {

  int conf, j, n, confnum, cadence, PART, confmax, fmc, other, end, num[3], skipnum, k, offset;
  char filename[100], line[300];
  double cm[3][3];
  float time4;
  float dumpf;
  FILE *in, *Cout, *Oout, *Neout, *chargefile;

  if(argc != 10) {
    fprintf(stderr,"Syntax Error: %s <input file> <charge file> <output file base> <# of configurations> <fm/c between ckpts> <max number of configuration steps> <# of particles> <# of configs to skip> <fm/c offset>\n",argv[0]);
    return 0;
  }

  chargefile=fopen(argv[2],"r");
  confnum=atoi(argv[4]);
  cadence=atoi(argv[5]);
  confmax=atoi(argv[6]);
  PART=atoi(argv[7]);
  skipnum=atoi(argv[8]);
  offset=atoi(argv[9]);

  sprintf(filename,"%s.C.out",argv[3]);
  Cout=fopen(filename,"w");
  sprintf(filename,"%s.O.out",argv[3]);
  Oout=fopen(filename,"w");
  sprintf(filename,"%s.Ne.out",argv[3]);
  Neout=fopen(filename,"w");

  float vel[confmax][PART][3], skip[PART][3];
  double pass, vaf[3][confmax], charge[PART];
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

  sprintf(filename,"%s",argv[1]);
  k=strlen(filename);
  open_xv4a_(filename,&PART,&time4,&skip[0][0],&skip[0][0],k);
  //  in=fopen(argv[1],"r");

  for(conf=0;conf<skipnum;conf++) {
    read_xv4a_(filename,&PART,&time4,&skip[0][0],&skip[0][0],k);
    /* fread(&dumpf,sizeof(float),1,in); */
    /* fread(&dumpf,sizeof(float),1,in); */
    /* fread(&dumpf,sizeof(float),1,in); */
    /* fread(&dumpf,sizeof(float),1,in); */
    /* for(j=0;j<PART;j++) */
    /*   fread(vel[conf%confmax][j],sizeof(float),3,in); */
    /* fread(&dumpf,sizeof(float),1,in); */
    printf("Skipping configuration %d\n",conf+1);
  }

  for(conf=0;conf<confnum;conf++) {
    printf("Configuration %d %e %e %e\n",conf+1,vaf[0][0],vaf[1][0],vaf[2][0]);
    if(load[conf]==0) {
      read_xv4a_(filename,&PART,&time4,&skip[0][0],&vel[conf%confmax][0][0],k);
      /* fread(&dumpf,sizeof(float),1,in); */
      /* fread(&dumpf,sizeof(float),1,in); */
      /* fread(&dumpf,sizeof(float),1,in); */
      /* fread(&dumpf,sizeof(float),1,in); */
      /* for(j=0;j<PART;j++) */
      /* 	fread(vel[conf%confmax][j],sizeof(float),3,in); */
      /* fread(&dumpf,sizeof(float),1,in); */
      load[conf]=1;
    }
    end=min(conf+confmax,confnum);
    for(other=conf;other<end;other++) {
      if(load[other]==0) {      
	read_xv4a_(filename,&PART,&time4,&skip[0][0],&vel[other%confmax][0][0],k);
	/* fread(&dumpf,sizeof(float),1,in); */
	/* fread(&dumpf,sizeof(float),1,in); */
	/* fread(&dumpf,sizeof(float),1,in); */
	/* fread(&dumpf,sizeof(float),1,in); */
	/* for(j=0;j<PART;j++) */
	/*   fread(vel[conf%confmax][j],sizeof(float),3,in); */
	/* fread(&dumpf,sizeof(float),1,in); */
	load[other]=1;
      }
      for(j=0;j<PART;j++) {
	pass=vel[conf%confmax][j][0]*vel[other%confmax][j][0]+vel[conf%confmax][j][1]*vel[other%confmax][j][1]+vel[conf%confmax][j][2]*vel[other%confmax][j][2];
	if(charge[j]<7.0)
	  vaf[0][other-conf]+=pass;
	else if(charge[j]<9.0)
	  vaf[1][other-conf]+=pass;
	else
	  vaf[2][other-conf]+=pass;
	/* 	if(charge[j]<7.0) */
	/* 	  vaf[0][other-conf]+=(vel[conf%confmax][j][0]-cm[0][0])*(vel[other%confmax][j][0]-cm[0][0])+(vel[conf%confmax][j][1]-cm[0][1])*(vel[other%confmax][j][1]-cm[0][1])+(vel[conf%confmax][j][2]-cm[0][2])*(vel[other%confmax][j][2]-cm[0][2]); */
	/* 	else if(charge[j]<9.0) */
	/* 	  vaf[1][other-conf]+=(vel[conf%confmax][j][0]-cm[1][0])*(vel[other%confmax][j][0]-cm[1][0])+(vel[conf%confmax][j][1]-cm[1][1])*(vel[other%confmax][j][1]-cm[1][1])+(vel[conf%confmax][j][2]-cm[1][2])*(vel[other%confmax][j][2]-cm[1][2]); */
	/* 	else  */
	/* 	  vaf[2][other-conf]+=(vel[conf%confmax][j][0]-cm[2][0])*(vel[other%confmax][j][0]-cm[2][0])+(vel[conf%confmax][j][1]-cm[2][1])*(vel[other%confmax][j][1]-cm[2][1])+(vel[conf%confmax][j][2]-cm[2][2])*(vel[other%confmax][j][2]-cm[2][2]); */
      } 
    }
  }
  //  fclose(in);
  close_xv4a_(filename,&PART,&time4,&skip[0][0],&skip[0][0],k);
  
  for(j=0;j<confmax;j++) {
    vaf[0][j]/=(num[0]);
    fprintf(Cout,"%d %e\n",j*cadence,vaf[0][j]);
    vaf[1][j]/=(num[1]);
    fprintf(Oout,"%d %e\n",j*cadence,vaf[1][j]);
    vaf[2][j]/=(num[2]);
    fprintf(Neout,"%d %e\n",j*cadence,vaf[2][j]);
  }
  
}
