#include <stdio.h>
#include <math.h>
#include <strings.h>

int main(int argc, char *argv[]) {

  int PART, slices, j, n, flag;
  double rho, BOXSIZE, sliced;
  char line[300], comp[1];
  FILE *in, *out;

  if(argc != 5) {
    fprintf(stderr,"Syntax Error: %s <input file> <output file> <# of slices> <density>\n",argv[0]);
    return 1;
  }

  in=fopen(argv[1],"r");
  out=fopen(argv[2],"w");
  slices=atoi(argv[3]);
  sscanf(argv[4],"%lf",&rho);

  fgets(line,300,in);
  sscanf(line,"%d",&PART);
  fgets(line,300,in);

  BOXSIZE=pow(PART/rho,1.0/3.0);
  sliced=BOXSIZE/slices;

  double pos[3][PART];
  int bins[slices][3];

  for(n=0;n<slices;n++)
    for(j=0;j<3;j++)
      bins[n][j]=0;

  for(j=0;j<PART;j++) {
    fgets(line,300,in);
    sscanf(line,"%s %lf %lf %lf\n",comp,&pos[0][j],&pos[1][j],&pos[2][j]);
    flag=0;
    for(n=0;n<slices;n++) {
      if(flag==0) {
	if(pos[2][j]<(n+1)*sliced) {
	  flag=1;
	  if(strcmp(comp,"C")==0)
	    bins[n][0]++;
	  else if(strcmp(comp,"O")==0)
	    bins[n][1]++;
	  else
	    bins[n][2]++;
	}
      }
    }
  }
  for(n=0;n<slices;n++)
    fprintf(out,"%e %d %d %d\n",n*sliced,bins[n][0],bins[n][1],bins[n][2]);
}
