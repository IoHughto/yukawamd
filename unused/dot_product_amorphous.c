#include <stdio.h>
#include <math.h>

#define PI 3.141592653589

int main(int argc, char *argv[]) {

  int PART, j, binnum, n, tflag, pflag, peaknum, peak, thresh, tPART, k;
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

  double pos[3][PART], theta[PART], phi[PART], max[14], r, mean, stdev, meanmax, e[3], trot[14], prot[14], x, y, z, c[3], cp, ct, s[3], sp, st, emax[3], de=0.1;
  for(j=0;j<PART;j++) {
    fgets(line,300,in);
    sscanf(line,"%lf %lf",&theta[j],&phi[j]);
  }
  tPART=PART/14;
  for(j=0;j<tPART;j++) {
    e[0]=0.0;
    e[1]=0.0;
    e[2]=0.0;
    meanmax=0.0;
    while(e[0]<2*PI) {
      e[1]=0.0;
      while(e[1]<PI) {
	  e[2]=0.0;
	  while(e[2]<2*PI) {
	      for(n=0;n<14;n++) {
		cp=cos(basis[1][n]);
		ct=cos(basis[0][n]);
		c[0]=cos(e[0]);
		c[1]=cos(e[1]);
		c[2]=cos(e[2]);
		sp=sin(basis[1][n]);
		st=sin(basis[0][n]);
		s[0]=sin(e[0]);
		s[1]=sin(e[1]);
		s[2]=sin(e[2]);
		x=cp*c[0]*c[1]*st-ct*s[0]+c[0]*sp*st*s[1];
		y=ct*c[0]*s[2]+cp*st*(c[1]*s[0]*s[2]-c[2]*s[1])+sp*st*(c[1]*c[2]+s[0]*s[1]*s[2]);
		z=ct*c[0]*c[2]+sp*st*(c[2]*s[0]*s[1]-c[1]*s[2])+cp*st*(c[1]*c[2]*s[0]+s[1]*s[2]);
		trot[n]=acos(z);
		prot[n]=acos(x/sqrt(x*x+y*y));
		if(y<0)
		  prot[n]*=-1.0;
	      }
	      for(k=0;k<14;k++) {
		max[k]=0.0;
		for(n=0;n<14;n++) {
		  r=cos(prot[n])*sin(trot[n])*cos(phi[14*j+k])*sin(theta[14*j+k])+sin(prot[n])*sin(trot[n])*sin(phi[14*j+k])*sin(theta[14*j+k])+cos(trot[n])*cos(theta[14*j+k]);
		  if(r>max[k])
		    max[k]=r;
		}
	      }
	      mean=0.0;
	      for(k=0;k<14;k++)
		mean+=max[k];
	      mean/=14;
	      if(mean>meanmax) {
		emax[0]=e[0];
		emax[1]=e[1];
		emax[2]=e[2];
		meanmax=mean;
	      }
	      e[2]+=de;
	    }
	  e[1]+=de;
	}
      e[0]+=de;
    }
    printf("%e %e %e\n",emax[0],emax[1],emax[2]);
  }
}

