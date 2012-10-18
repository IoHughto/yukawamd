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

  if(argc != 6) {
    fprintf(stderr,"Syntax Error: %s <solid file> <liquid file> <output config> <# of particles> <density in 1/fm^3>\n",argv[0]);
    return 0;
  }
 
  PART=atoi(argv[4]);
  sscanf(argv[5],"%lf",&rho);

  newPART=8*PART;

  BOXSIZE=pow(PART/rho,1.0/3.0);

  double pos[PART][3], vel[PART][3];
  double newpos[8*PART][3], newvel[8*PART][3];

  printf("here 1\n");

  // Solid bottom
  sprintf(filename,"%s",argv[1]);
  k=strlen(filename);
  read_xv8_(filename,&PART,&time8,&pos[0][0],&vel[0][0],k);

  printf("here 2\n");

  for(j=0;j<PART;j++)
    for(x=0;x<2;x++)
      for(y=0;y<2;y++) {
	newpos[x*PART+2*y*PART+j][0]=pos[j][0]+x*BOXSIZE;
	newpos[x*PART+2*y*PART+j][1]=pos[j][1]+y*BOXSIZE;
	newpos[x*PART+2*y*PART+j][2]=pos[j][2];
	newvel[x*PART+2*y*PART+j][0]=vel[j][0];
	newvel[x*PART+2*y*PART+j][1]=vel[j][1];
	newvel[x*PART+2*y*PART+j][2]=vel[j][2];
      }

  printf("here 3\n");

  // Liquid top
  sprintf(filename,"%s",argv[2]);
  k=strlen(filename);
  read_xv8_(filename,&PART,&time8,&pos[0][0],&vel[0][0],k);

  printf("here 4\n");

  for(j=0;j<PART;j++)
    for(x=0;x<2;x++)
      for(y=0;y<2;y++) {
	newpos[x*PART+2*y*PART+j+4*PART][0]=pos[j][0]+x*BOXSIZE;
	newpos[x*PART+2*y*PART+j+4*PART][1]=pos[j][1]+y*BOXSIZE;
	newpos[x*PART+2*y*PART+j+4*PART][2]=pos[j][2]+BOXSIZE;
	newvel[x*PART+2*y*PART+j+4*PART][0]=vel[j][0];
	newvel[x*PART+2*y*PART+j+4*PART][1]=vel[j][1];
	newvel[x*PART+2*y*PART+j+4*PART][2]=vel[j][2];
      }

  printf("here 5\n");

  time8=0.0;
  sprintf(filename,"%s",argv[3]);
  k=strlen(filename);
  write_xv8_(filename,&newPART,&time8,&newpos[0][0],&newvel[0][0],k);

  return 0;

}
