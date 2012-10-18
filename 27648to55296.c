#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

void read_xv8_(char*,int*,double*,double*,double*,int);
void write_xv8_(char*,int*,double*,double*,double*,int);

int main(int argc, char *argv[]) {

  double BOXSIZE, rho, time8;
  int PART, j, x, y, newPART, k;
  char filename[100], line[300];
  FILE *solid, *liquid, *out;

  if(argc != 5) {
    fprintf(stderr,"Syntax Error: %s <solid file> <liquid file> <output config> <density in 1/fm^3>\n",argv[0]);
    return 0;
  }
 
  PART=27648;
  sscanf(argv[4],"%lf",&rho);

  newPART=2*PART;

  BOXSIZE=pow(PART/rho,1.0/3.0);

  double pos[PART][3], vel[PART][3];
  double newpos[2*PART][3], newvel[2*PART][3];

  // Solid bottom
  sprintf(filename,"%s",argv[1]);
  k=strlen(filename);
  read_xv8_(filename,&PART,&time8,&pos[0][0],&vel[0][0],k);

  for(j=0;j<PART;j++) {
        newpos[j][0]=pos[j][0];
        newpos[j][1]=pos[j][1];
        newpos[j][2]=pos[j][2];
        newvel[j][0]=vel[j][0];
        newvel[j][1]=vel[j][1];
        newvel[j][2]=vel[j][2];
  }

  // Liquid top
  sprintf(filename,"%s",argv[2]);
  k=strlen(filename);
  read_xv8_(filename,&PART,&time8,&pos[0][0],&vel[0][0],k);

  for(j=0;j<PART;j++) {
	newpos[PART+j][0]=pos[j][0];
	newpos[PART+j][1]=pos[j][1];
	newpos[PART+j][2]=pos[j][2]+BOXSIZE;
	newvel[PART+j][0]=vel[j][0];
	newvel[PART+j][1]=vel[j][1];
	newvel[PART+j][2]=vel[j][2];
 }

  time8=0.0;
  sprintf(filename,"%s",argv[3]);
  k=strlen(filename);
  write_xv8_(filename,&newPART,&time8,&newpos[0][0],&newvel[0][0],k);

  return 0;

}
