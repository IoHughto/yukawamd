#include <stdio.h>
#include <math.h>

#define PI 3.141592653589

int main(int argc, char *argv[]) {

  int PART, j, binnum, n, tflag, pflag, peaknum, peak, thresh;
  double bind;
  char line[300], dumpc, filename[100];
  FILE *in, *out[2], *basisin;

  if(argc!=6) {
    fprintf(stderr,"Syntax Error: %s <input (theta/phi) file> <basis vector file> <output file> <# of bins> <# of particles>\n",argv[0]);
    return 1;
  }

  in=fopen(argv[1],"r");
  basisin=fopen(argv[2],"r");
  out[0]=fopen(argv[3],"w");
  sprintf(filename,"%s.full",argv[3]);
  out[1]=fopen(filename,"w");
  binnum=atoi(argv[4]);
  PART=atoi(argv[5]);

  bind=PI/binnum;

  int bins[binnum];
  double basis[2][14];

  for(n=0;n<binnum;n++) {
    bins[n]=0;
  }

  for(n=0;n<14;n++) {
    fgets(line,300,basisin);
    sscanf(line,"%lf %lf",&basis[0][n],&basis[1][n]);
  }

  double pos[3][PART], theta[PART], phi[PART], max[14], r, mean, stdev;
  for(j=0;j<PART;j++) {
    fgets(line,300,in);
    sscanf(line,"%lf %lf",&theta[j],&phi[j]);
    max[j%14]=0.0;
    for(n=0;n<14;n++) {
      r=cos(basis[1][n])*sin(basis[0][n])*cos(phi[j])*sin(theta[j])+sin(basis[1][n])*sin(basis[0][n])*sin(phi[j])*sin(theta[j])+cos(basis[0][n])*cos(theta[j]);
      if(r>max[j%14])
	max[j%14]=r;
    }
    if((j+1)%14==0) {
      for(n=0;n<14;n++) 
       	fprintf(out[1],"%e ",max[n]); 
      fprintf(out[1],"\n"); 
      mean=0.0;
      for(n=0;n<14;n++) 
	mean+=max[n];
      mean/=14.0;
      stdev=0.0;
      for(n=0;n<14;n++)
	stdev+=(mean-max[n])*(mean-max[n]);
      stdev/=13.0;
      stdev=sqrt(stdev);
      if(stdev<0.024) 
       	fprintf(out[0],"C\n"); 
      else {
	if(mean>0.955)
	  fprintf(out[0],"C\n");
	else
	  fprintf(out[0],"O\n");
      } 
      /* if(mean>0.955) {  */
      /*  	if(stdev<0.025)  */
      /*  	  fprintf(out[0],"C\n");  */
      /*  	else  */
      /*  	  fprintf(out[0],"N\n");  */
      /* }  */
      /* else {  */
      /*  	if(stdev>0.025)  */
      /*  	  fprintf(out[0],"O\n");  */
      /*  	else  */
      /*  	  fprintf(out[0],"F\n");  */
      /* }  */
    }
  }
}
