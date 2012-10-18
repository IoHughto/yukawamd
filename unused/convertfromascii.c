#include <stdio.h>

int main(int argc, char *argv[]) {

  int PART, j;
  char line[500];
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
    fgets(line,500,in);
    sscanf(line,"%lf %lf %lf %lf %lf %lf\n",&pos[j][0],&pos[j][1],&pos[j][2],&vel[j][0],&vel[j][1],&vel[j][2]);
    fwrite(pos[j],sizeof(double),sizeof(pos[j])/sizeof(double),out);
    fwrite(vel[j],sizeof(double),sizeof(vel[j])/sizeof(double),out);
  }
}
