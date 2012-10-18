#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {

  int j, n, PART;
  double time=0.0;
  float dump=0.0;
  FILE *in, *out;

  if(argc != 4) {
    fprintf(stderr,"Syntax Error: %s <input> <output> <# of particles>\n",argv[0]);
    return 0;
  }

  in=fopen(argv[1],"r");
  out=fopen(argv[2],"w");

  PART=atoi(argv[3]);

  double pos[3*PART], vel[3*PART];

  fwrite(&dump,sizeof(float),1,out);
  fwrite(&time,sizeof(double),1,out);
  fwrite(&dump,sizeof(float),1,out);
  fwrite(&dump,sizeof(float),1,out);

  fread(pos,sizeof(double),3*PART,in);
  fread(vel,sizeof(double),3*PART,in);

  for(j=0;j<PART;j++) { 
    for(n=0;n<3;n++) 
      fwrite(&pos[3*j+n],sizeof(double),1,out);
    for(n=0;n<3;n++) 
      fwrite(&vel[3*j+n],sizeof(double),1,out);
  }
  fwrite(&dump,sizeof(float),1,out);

  fclose(in);
  fclose(out);
}
