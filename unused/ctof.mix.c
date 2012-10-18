#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {

  unsigned char a;
  int j, n, PART;
  double time=0.0, trans;
  float dump=0.0;
  FILE *in[2], *out;

  if(argc != 5) {
    fprintf(stderr,"Syntax Error: %s <input 1> <input 2> <output> <# of particles>\n",argv[0]);
    return 0;
  }

  in[0]=fopen(argv[1],"r");
  in[1]=fopen(argv[2],"r");
  out=fopen(argv[3],"w");

  PART=atoi(argv[4]);

  double pos[3*PART], vel[3*PART];

  for(j=0;j<20;j++) {
    fread(&a,sizeof(char),1,in[0]);
    fwrite(&a,sizeof(char),1,out);
  }

  fread(pos,sizeof(double),3*PART,in[1]);
  fread(vel,sizeof(double),3*PART,in[1]);

  for(j=0;j<6*3456;j++)
    fread(&trans,sizeof(double),1,in[0]);

  for(j=0;j<PART;j++) { 
    for(n=0;n<3;n++) 
      fwrite(&pos[3*j+n],sizeof(double),1,out);
    for(n=0;n<3;n++) 
      fwrite(&vel[3*j+n],sizeof(double),1,out);
  }
  for(j=0;j<4;j++) {
    fread(&a,sizeof(char),1,in[0]);
    fwrite(&a,sizeof(char),1,out);
  }

  fclose(in[0]);
  fclose(in[1]);
  fclose(out);
}
