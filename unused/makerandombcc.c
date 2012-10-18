#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.141592653589

double gaussran(double *pvx, double *pvy, double *pvz, double T, double m);

int main(int argc, char *argv[])
{
  srand(time(NULL)<<8);

  int j, k, l, PART, n, dump, num, m;
  double rho, numd, T, theta, phi, psi, x, y, z, r, BOXSIZE, temp[3];
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

  double pos[2*PART][3],vel[2*PART][3], charge[2*PART], mass[2*PART];

  for(j=0;j<PART;j++) {
    fgets(line,300,chargefile);
    sscanf(line,"%d %lf %lf\n",&dump,&charge[j],&mass[j]);
    gaussran(&vel[j][0],&vel[j][1],&vel[j][2],T,mass[j]);
  }

  numd=pow((double)PART/2.0,1.0/3.0)+0.1;
 
  BOXSIZE=pow(PART/rho,1.0/3.0);

  num=numd;

  stretch=pow(2.0/rho,1.0/3.0);

  if(num*num*num*2 != PART) {
    fprintf(stderr,"%d particles do not make a complete BCC system.\n",PART);
  }

  theta=(double)rand()/RAND_MAX*PI/2.0;
  phi=(double)rand()/RAND_MAX*PI/2.0;
  psi=(double)rand()/RAND_MAX*PI/2.0;
  x=(double)rand()/RAND_MAX*stretch;
  y=(double)rand()/RAND_MAX*stretch;
  z=(double)rand()/RAND_MAX*stretch;
  
  printf("%lf %lf %lf %lf %lf %lf\n",theta,phi,psi,x,y,z);

  n=0;
  for(j=-2*num;j<2*num;j++) {
    for(k=-2*num;k<2*num;k++) {
      for(l=-2*num;l<2*num;l++) {
	temp[0]=stretch*(double)j;
	temp[1]=stretch*(double)k;
	temp[2]=stretch*(double)l;
	pos[n][0]=temp[0]*cos(theta)*cos(psi)+temp[1]*(sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi))+temp[2]*(sin(phi)*sin(psi)+cos(phi)*sin(theta)*cos(psi));
	pos[n][1]=temp[0]*cos(theta)*sin(psi)+temp[1]*(sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi))+temp[2]*(-1.0*sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi));
	pos[n][2]=-1.0*temp[0]*sin(theta)+temp[1]*sin(phi)*cos(theta)+temp[2]*cos(phi)*cos(theta);
	if(pos[n][0]>0.0)
	  if(pos[n][0]<BOXSIZE)
	    if(pos[n][1]>0.0) 
	      if(pos[n][1]<BOXSIZE)
		if(pos[n][2]>0.0)
		  if(pos[n][2]<BOXSIZE)
		    n++;
      }
    }
  }  

  for(j=-2*num;j<2*num;j++) {
    for(k=-2*num;k<2*num;k++) {
      for(l=-2*num;l<2*num;l++) {
	temp[0]=stretch*(double)j+stretch/2.0;
	temp[1]=stretch*(double)k+stretch/2.0;
	temp[2]=stretch*(double)l+stretch/2.0;
	pos[n][0]=temp[0]*cos(theta)*cos(psi)+temp[1]*(sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi))+temp[2]*(sin(phi)*sin(psi)+cos(phi)*sin(theta)*cos(psi));
	pos[n][1]=temp[0]*cos(theta)*sin(psi)+temp[1]*(sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi))+temp[2]*(-1.0*sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi));
	pos[n][2]=-1.0*temp[0]*sin(theta)+temp[1]*sin(phi)*cos(theta)+temp[2]*cos(phi)*cos(theta);
	if(pos[n][0]>0.0)
	  if(pos[n][0]<BOXSIZE)
	    if(pos[n][1]>0.0) 
	      if(pos[n][1]<BOXSIZE)
		if(pos[n][2]>0.0)
		  if(pos[n][2]<BOXSIZE)
		    n++;
      }
    }
  }

  printf("%d %d\n",PART,n);

  fprintf(out,"%10.8lf\n",0.0);

  for(j=0;j<n;j++) 
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
