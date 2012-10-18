#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char *argv[]) {
  srand(time(NULL)<<8);

  int PART, j;
  double value, charge, OCP, TCP;
  FILE *chargefile;

  if(argc != 4) {
    fprintf(stderr,"Syntax Error: %s <# of particles> <X_1> <X_2>\n",argv[0]);
    return 0;
  }

  chargefile=fopen("mass.config","w");
  PART=atoi(argv[1]);
  sscanf(argv[2],"%lf",&OCP);
  sscanf(argv[3],"%lf",&TCP);

  for(j=0;j<PART;j++) {
    value=(double)rand()/RAND_MAX;

    if(value<OCP)
      charge=10.0;
    else if(value<OCP+TCP)
      charge=6.0;
    else
      charge=8.0;

    fprintf(chargefile,"%lf\n",charge);
  }

  return 0;
}
