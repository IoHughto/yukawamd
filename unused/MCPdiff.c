#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

void read_xv8_(char*,int*,double*,double*,double*,int);
void write_xv8_(char*,int*,double*,double*,double*,int);

int main(int argc, char *argv[]) {
  int j, k, l, n, PART, qflag, flag, difference, newnum;
  double rcut, rho, BOXSIZE, time8, diff[3], r, delta, newZ, newmass, sfac, impfrac;
  FILE *inion, *oution;
  char filename[100], line[300];

  if(argc != 9) {
    fprintf(stderr,"Syntax Error: %s <ckpt in> <ion.dat in> <density> <new Z> <fraction of impurity> <ckpt out> <ion.dat out> <Randomize positions?>\n",argv[0]);
    return 1;
  }
  srand ( time(NULL) );

  inion=fopen(argv[2],"r");
  sscanf(argv[3],"%lf",&rho);
  sscanf(argv[4],"%lf",&newZ);
  sscanf(argv[5],"%lf",&impfrac);
  qflag=atoi(argv[8]);

  fgets(line,300,inion);
  sscanf(line,"%d",&PART);
  
  newnum=(int)(PART*impfrac);

  int cutparticles[PART], newparticles[PART], ncut, nnew, newion[newnum];
  double pos[PART][3], vel[PART][3], num;
  double crystal[PART][3];
  double charge[PART], mass[PART];

  BOXSIZE=pow(PART/rho,1.0/3.0);
  delta=pow(2.0/rho,1.0/3.0);
  num=pow(PART/2.0,1.0/3.0);

  j=0;

  while(fgets(line,300,inion)!=NULL) {
    sscanf(line,"%d %lf %lf\n",&k,&charge[j],&mass[j]);
    j++;
  }

  switch((int)newZ) {
  case 2:
    newmass=4.0;
    break;
  case 4:
    newmass=8.0;
    break;
  case 6:
    newmass=12.0;
    break;
  case 8:
    newmass=16.0;
    break;
  case 10:
    newmass=22.0;
    break;
  case 12:
    newmass=24.0;
    break;
  case 14:
    newmass=28.0;
    break;
  case 16:
    newmass=32.0;
    break;
  case 18:
    newmass=40.0;
    break;
  case 20:
    newmass=40.0;
    break;
  case 22:
    newmass=48.0;
    break;
  case 24:
    newmass=52.0;
    break;
  case 26:
    newmass=56.0;
    break;
  case 28:
    newmass=58.0;
    break;
  case 30:
    newmass=66.0;
    break;
  case 32:
    newmass=72.0;
    break;
  case 34:
    newmass=80.0;
    break;
  case 36:
    newmass=84.0;
    break;
  case 38:
    newmass=88.0;
    break;
  case 40:
    newmass=92.0;
    break;
  default:
    fprintf(stderr,"Ion error: New Z out of bounds.\n");
    return 1;
    break;
  }

  sprintf(filename,"%s",argv[1]);
  k=strlen(filename);
  read_xv8_(filename,&PART,&time8,&pos[0][0],&vel[0][0],k);

  for(j=0;j<newnum;j++) {
    newion[j]=rand()%PART;
    sfac=sqrt(mass[newion[j]]/newmass);
    mass[newion[j]]=newmass;
    charge[newion[j]]=newZ;
    for(n=0;n<3;n++)
      vel[newion[j]][n]*=sfac;
    if(qflag==1) {
      for(n=0;n<3;n++) 
	pos[newion[j]][n]=rand()/RAND_MAX;
    }
  }

  sprintf(filename,"%s",argv[6]);
  k=strlen(filename);
  write_xv8_(filename,&PART,&time8,&pos[0][0],&vel[0][0],k);
  
  oution=fopen(argv[7],"w");
  fprintf(oution,"%d\n",PART);
  for(j=0;j<PART;j++) {
    fprintf(oution,"%d %lf %lf\n",j,charge[j],mass[j]);
  }
}
