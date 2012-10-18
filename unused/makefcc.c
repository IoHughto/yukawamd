#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.141592653589

double gaussran(double *pvx, double *pvy, double *pvz, double T, double m);

int main(int argc, char *argv[])
{
  int j, k, l, PART, n, dump, num;
  double rho, numd, T;
  double stretch;
  char line[300];
  FILE *chargefile, *out;

  if(argc != 5) {
    fprintf(stderr,"Syntax Error: %s <output file> <charge.dat file> <density in 1/fm^3> <initial T>\n",argv[0]);
    return 0;
  }

  out=fopen(argv[1],"w");
  chargefile=fopen(argv[2],"r");
  sscanf(argv[3],"%lf",&rho);
  sscanf(argv[4],"%lf",&T);

  fgets(line,300,chargefile);
  sscanf(line,"%d",&PART);

  double pos[PART][3],vel[PART][3], charge[PART], mass[PART];

  for(j=0;j<PART;j++) {
    fgets(line,300,chargefile);
    sscanf(line,"%d %lf %lf\n",&dump,&charge[j],&mass[j]);
    mass[j]*=931.0;
    gaussran(&vel[j][0],&vel[j][1],&vel[j][2],T,mass[j]);
  }

  numd=pow((double)PART/4.0,1.0/3.0)+0.1;

  num=numd;

  stretch=pow(4.0/rho,1.0/3.0);

  if(num*num*num*4 != PART) {
    fprintf(stderr,"%d particles do not make a complete BCC system.\n",PART);
    return 2;
  }

  n=0;
  for(j=0;j<num;j++) {
    for(k=0;k<num;k++) {
      for(l=0;l<num;l++) {
	pos[n][0]=stretch*(double)j;
	pos[n][1]=stretch*(double)k;
	pos[n][2]=stretch*(double)l;
	n++;
      }
    }
  }  
  for(j=0;j<num;j++) {
    for(k=0;k<num;k++) {
      for(l=0;l<num;l++) {
	pos[n][0]=stretch*(double)j;
	pos[n][1]=stretch*(double)k+stretch/2.0;
	pos[n][2]=stretch*(double)l+stretch/2.0;
	n++;
      }
    }
  }  
  for(j=0;j<num;j++) {
    for(k=0;k<num;k++) {
      for(l=0;l<num;l++) {
	pos[n][0]=stretch*(double)j+stretch/2.0;
	pos[n][1]=stretch*(double)k+stretch/2.0;
	pos[n][2]=stretch*(double)l;
	n++;
      }
    }
  }  
  for(j=0;j<num;j++) {
    for(k=0;k<num;k++) {
      for(l=0;l<num;l++) {
	pos[n][0]=stretch*(double)j+stretch/2.0;
	pos[n][1]=stretch*(double)k;
	pos[n][2]=stretch*(double)l+stretch/2.0;
	n++;
      }
    }
  }

  fprintf(out,"%10.8lf\n",0.0);

  for(j=0;j<PART;j++) 
    fprintf(out,"%16.14e %16.14e %16.14e %16.14e %16.14e %16.14e\n",pos[j][0],pos[j][1],pos[j][2],vel[j][0],vel[j][1],vel[j][2]);

  fclose(out);

}

double gaussran(double *pvx, double *pvy, double *pvz, double T, double m) {
  int n;
  double x2, y2, x, y, z, theta, phi;
  
  // Rejection method for sampling under a Boltzmann dist.
  do{
    x2=(double)rand()/RAND_MAX;
    y2=4.0/exp(1.0)*sqrt(m/(2.0*PI*T))*(double)rand()/(RAND_MAX);
  } while(y2>(4*PI*(m/(2.0*PI*T))*sqrt(m/(2.0*PI*T))*(double)x2*(double)x2*exp(-m*(double)x2*(double)x2/(2.0*T))));
  
  // Choose phi between 0 and 2pi, and cos(theta) between -1 and 1
  theta=2*(double)rand()/RAND_MAX-1;
  phi=2*PI*(double)rand()/RAND_MAX;
  
  // Give Cartesian velocities back
  *pvx=x2*cos(phi)*sqrt(1-theta*theta);
  *pvy=x2*sin(phi)*sqrt(1-theta*theta);
  *pvz=x2*theta;
  
  return 0;
}
