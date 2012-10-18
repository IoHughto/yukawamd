#include <stdio.h>

int main(int argc, char *argv[]) {
  int n, j, PART;
  FILE *in, *out;

  if(argc != 4) {
    fprintf(stderr,"Syntax Error: %s <input file> <output file> <# of particles>\n",argv[0]);
    return 1;
  }
  
  in=fopen(argv[1],"r");
  out=fopen(argv[2],"w");

  PART=atoi(argv[3]);

  fprintf(out,"%15.3f\n",0.0);
  
  double pos[3*PART], vel[3*PART];

  fread(pos,sizeof(double),3*PART,in);
  fread(vel,sizeof(double),3*PART,in);

  for(j=0;j<PART;j++) {
    for(n=0;n<3;n++) {
      fprintf(out," % 18.16e",pos[3*j+n]);
    }
    for(n=0;n<3;n++) {
      fprintf(out," % 18.16e",vel[3*j+n]);
    }
    fprintf(out,"\n");
  }
  fclose(in);
  fclose(out);
}
