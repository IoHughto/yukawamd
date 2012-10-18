#include <stdio.h>
#include <math.h>

int main(int argc, char *argv[]) {
  int PART, j, k, l;
  char line[300], dumpc;
  FILE *in, *out;

  if(argc != 3) {
    fprintf(stderr,"Syntax Error: %s <input file> <output file>\n",argv[0]);
    return 1;
  }

  in=fopen(argv[1],"r");
  out=fopen(argv[2],"w");

  fgets(line,300,in);
  sscanf(line,"%d",&PART);
  fgets(line,300,in);

  PART/=14;

  double pos[3][14], theta, r[14];

  for(j=0;j<PART;j++) {
    for(k=0;k<14;k++) {
      fgets(line,300,in);
      sscanf(line,"%c %lf %lf %lf\n",&dumpc,&pos[0][k],&pos[1][k],&pos[2][k]);
      r[k]=sqrt(pos[0][k]*pos[0][k]+pos[1][k]*pos[1][k]+pos[2][k]*pos[2][k]);
    }
    for(k=0;k<14;k++)
      for(l=0;l<k;l++) {
	theta=(pos[0][k]*pos[0][l]+pos[1][k]*pos[1][l]+pos[2][k]*pos[2][l])/(r[k]*r[l]);
	fprintf(out,"%lf ",theta);
      }
    fprintf(out,"\n");
  }
}
