#include <stdio.h>

int main(int argc, char *argv[]) {
 
  int j, n, fmc, offset, cadence, confnum;
  char line[300], filename[100];
  FILE *in, *out;

  if(argc != 6) {
    fprintf(stderr,"Syntax Error: %s <file base> <output file> <fm/c cadence> <fm/c offset> <# of configs>\n",argv[0]);
    return 0;
  }

  out=fopen(argv[2],"w");
  cadence=atoi(argv[3]);
  offset=atoi(argv[4]);
  confnum=atoi(argv[5]);

  double dumpd[3], vel[3];

  for(n=0;n<confnum;n++) {
    fmc=(n+1)*cadence+offset;
    sprintf(filename,"%s.%011d",argv[1],fmc);
    in=fopen(filename,"r");
    fgets(line,300,in);
    while(fgets(line,300,in) != NULL) {
      sscanf(line,"%lf %lf %lf %lf %lf %lf",&dumpd[0],&dumpd[1],&dumpd[2],&vel[0],&vel[1],&vel[2]);
      fwrite(vel,3,sizeof(double),out);
    }
    fclose(in);
  }
}
