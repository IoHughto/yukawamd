#include <stdio.h>

int main(int argc, char *argv[]) {

  int filenum, total, PART, n, j, perfile;
  char filename[50];
  FILE *in, *out;

  if(argc != 5) {
    fprintf(stderr,"Syntax Error: %s <input file> <# of files> <# of configs> <# of particles>\n",argv[0]);
    return 1;
  }

  in=fopen(argv[1],"r");

  filenum=atoi(argv[2]);
  total=atoi(argv[3]);
  PART=atoi(argv[4]);

  double pos[3*PART], vel[3*PART];

  perfile=total/filenum;

  for(n=0;n<filenum;n++) {
    sprintf(filename,"%d.out",n);
    out=fopen(filename,"w");
    for(j=0;j<perfile;j++){
      fread(pos,sizeof(double),3*PART,in);
      fread(vel,sizeof(double),3*PART,in);
      fwrite(pos,sizeof(double),3*PART,out);
      fwrite(vel,sizeof(double),3*PART,out);
    }
    fclose(out);
  }

  return 0;
}
