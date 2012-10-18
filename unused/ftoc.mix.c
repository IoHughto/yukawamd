#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {

  int j, n, PART;
  double time=0.0, trans;
  float dump=0.0;
  FILE *in, *out;

  if(argc != 4) {
    fprintf(stderr,"Syntax Error: %s <input> <output> <# of particles>\n",argv[0]);
    return 0;
  }

  in=fopen(argv[1],"r");
  out=fopen(argv[2],"w");

  PART=atoi(argv[3]);

  double pos[PART][3], vel[PART][3];
  
  fread(&dump,sizeof(float),1,in);
  fread(&time,sizeof(double),1,in);
  fread(&dump,sizeof(float),1,in);
  fread(&dump,sizeof(float),1,in);

  for(j=0;j<PART;j++) {
    fread(pos[j],sizeof(double),3,in);
    fread(vel[j],sizeof(double),3,in);
  }
  for(j=0;j<PART;j++) 
    fwrite(pos[j],sizeof(double),3,out);
  for(j=0;j<PART;j++) 
    fwrite(vel[j],sizeof(double),3,out);
  
  fread(&dump,sizeof(float),1,in);

  fclose(in);
  fclose(out);
}
