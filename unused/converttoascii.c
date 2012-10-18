#include <stdio.h>

int main(int argc, char *argv[]) {

  int PART, j;
  FILE *in, *out;

  if(argc != 4) {
    fprintf(stderr,"Syntax Error: %s <# of particles> <input file> <output file>\n",argv[0]);
    return 0;
  }

  PART=atoi(argv[1]);
  in=fopen(argv[2],"r");
  out=fopen(argv[3],"w");

  double pos[PART][3], vel[PART][3];

  for(j=0;j<PART;j++) {
    fread(pos[j],sizeof(double),sizeof(pos[j])/sizeof(double),in);
    fread(vel[j],sizeof(double),sizeof(vel[j])/sizeof(double),in);
    fprintf(out,"%32.30lf %32.30lf %32.30lf %32.30lf %32.30lf %32.30lf\n",pos[j][0],pos[j][1],pos[j][2],vel[j][0],vel[j][1],vel[j][2]);
  }
}
