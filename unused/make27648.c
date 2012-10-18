#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Define constants
#define PI 3.141592653589
#define ALPHA 0.00729735
#define HBARC 197.32693

double gaussran(double *pvx, double *pvy, double *pvz, double T, double m);
double ran(double *pvx, double *pvy, double *pvz, double BOXSIZE);

int main( int argc, char *argv[]) {
  int PART, j, n, dump, copy, x[3];
  double rho, BOXSIZE, T, dumpf;
  char line[300];
  FILE *par, *in, *out, *chin, *chout; 

  if(argc != 5) {
      fprintf(stderr,"Syntax Error: %s <input file> <output file> <charge input> <charge output>(you also need a parameters.in file)\n",argv[0]);
      return 0;
    }

  // Open up files for reading and writing
  in=fopen(argv[1],"r");
  par=fopen("parameters.in","r");
  out=fopen(argv[2],"w");
  chin=fopen(argv[3],"r");
  chout=fopen(argv[4],"w");

  fgets(line,300,par);
  fgets(line,300,par);
  sscanf(line,"%lf %d %lf",&rho,&PART,&dumpf);
  fgets(line,300,par);
  fgets(line,300,par);
  sscanf(line,"%lf %lf %lf",&dumpf,&dumpf,&dumpf);
  fgets(line,300,par);
  fgets(line,300,par);
  sscanf(line,"%d %lf %lf",&dump,&dumpf,&T);

  BOXSIZE=pow(PART/rho,1.0/3.0);

  double velocity[3*PART], position[3*PART], charge[PART], vel[8][3*PART], pos[8][3*PART];
  fread(position,sizeof(double),3*PART,in);
  fread(velocity,sizeof(double),3*PART,in);
  
  for(j=0;j<PART;j++) {
    fgets(line,300,chin);
    sscanf(line,"%lf",&charge[j]);
  }

  for(x[0]=0;x[0]<2;x[0]++) {
    for(x[1]=0;x[1]<2;x[1]++) {
      for(x[2]=0;x[2]<2;x[2]++) {
	if(x[2]<1) {
	  for(j=0;j<PART;j++) {
	    for(n=0;n<3;n++) {
	      pos[4*x[0]+2*x[1]+x[2]][3*j+n]=position[3*j+n]+x[n]*BOXSIZE;
	      vel[4*x[0]+2*x[1]+x[2]][3*j+n]=velocity[3*j+n];
	    }
	  }
	}
	else {
	  for(j=0;j<PART;j++) {
	    gaussran(&vel[4*x[0]+2*x[1]+x[2]][3*j],&vel[4*x[0]+2*x[1]+x[2]][3*j+1],&vel[4*x[0]+2*x[1]+x[2]][3*j+2],T,charge[j]);
	    ran(&pos[4*x[0]+2*x[1]+x[2]][3*j],&pos[4*x[0]+2*x[1]+x[2]][3*j+1],&pos[4*x[0]+2*x[1]+x[2]][3*j+2],BOXSIZE);
	    for(n=0;n<3;n++) {
	      pos[4*x[0]+2*x[1]+x[2]][3*j+n]=pos[4*x[0]+2*x[1]+x[2]][3*j+n]+x[n]*BOXSIZE;
	    }
	  }
	}
	printf("%lf %lf %lf\n",pos[4*x[0]+2*x[1]+x[2]][999],pos[4*x[0]+2*x[1]+x[2]][1000],pos[4*x[0]+2*x[1]+x[2]][1001]);
	printf("%d %d %d\n",x[0],x[1],x[2]);
	for(j=0;j<PART;j++)
	  fprintf(chout,"%lf\n",charge[j]);
      }
    }
  }
  fwrite(pos[0],sizeof(double),3*PART,out);
  fwrite(pos[1],sizeof(double),3*PART,out);
  fwrite(pos[2],sizeof(double),3*PART,out);
  fwrite(pos[3],sizeof(double),3*PART,out);
  fwrite(pos[4],sizeof(double),3*PART,out);
  fwrite(pos[5],sizeof(double),3*PART,out);
  fwrite(pos[6],sizeof(double),3*PART,out);
  fwrite(pos[7],sizeof(double),3*PART,out);
  fwrite(vel[0],sizeof(double),3*PART,out);
  fwrite(vel[1],sizeof(double),3*PART,out);
  fwrite(vel[2],sizeof(double),3*PART,out);
  fwrite(vel[3],sizeof(double),3*PART,out);
  fwrite(vel[4],sizeof(double),3*PART,out);
  fwrite(vel[5],sizeof(double),3*PART,out);
  fwrite(vel[6],sizeof(double),3*PART,out);
  fwrite(vel[7],sizeof(double),3*PART,out);
}

double gaussran(double *pvx, double *pvy, double *pvz, double T, double m) {
  int n;
  double x2, y2, x, y, z, theta, phi, mass;
  
  mass=2*931.494088*m;

  // Rejection method for sampling under a Boltzmann dist.
  do{
    x2=(double)rand()/RAND_MAX;
    y2=4.0/exp(1.0)*sqrt(mass/(2.0*PI*T))*(double)rand()/(RAND_MAX);
  } while(y2>(4*PI*(mass/(2.0*PI*T))*sqrt(mass/(2.0*PI*T))*(double)x2*(double)x2*exp(-mass*(double)x2*(double)x2/(2.0*T))));
  
  // Choose phi between 0 and 2pi, and cos(theta) between -1 and 1
  theta=2*(double)rand()/RAND_MAX-1;
  phi=2*PI*(double)rand()/RAND_MAX;
  
  // Give Cartesian velocities back
  *pvx=x2*cos(phi)*sqrt(1-theta*theta);
  *pvy=x2*sin(phi)*sqrt(1-theta*theta);
  *pvz=x2*theta;
  
  return 0;
}

double ran(double *pvx, double *pvy, double *pvz, double BOXSIZE) {
  double w;


  // Put the particles in a box
  w=BOXSIZE*(double)rand()/RAND_MAX;
  *pvx=w;

  w=BOXSIZE*(double)rand()/RAND_MAX;
  *pvy=w;

  w=BOXSIZE*(double)rand()/RAND_MAX;
  *pvz=w;

  return 0;
}
