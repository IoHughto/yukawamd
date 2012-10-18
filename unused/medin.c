#include <stdio.h>
#include <math.h>

double f_l(double g);
double f_s(double g);
double f_lTCP(double g, double x1, double x2, double Z1, double Z2);
double f_sTCP(double g, double x1, double x2, double Z1, double Z2);
double Delta_g(double x, double rz);

int main(int argc, char *argv[]) {
  
  double Z[2], a[2], b[2], x[2], xmesh, norm, g, fTCP[3], A, min, gstart, gstop, gmesh;
  int cbin, gbin, i, j, k, phase, iter, confnum;
  
  if(argc != 5) {
    fprintf(stderr,"Syntax Error: %s <Gamma start> <Gamma stop> <Gamma mesh> <xmesh>\n",argv[0]);
  }

  sscanf(argv[1],"%lf",&gstart);
  sscanf(argv[2],"%lf",&gstop);
  sscanf(argv[3],"%lf",&gmesh);
  sscanf(argv[4],"%lf",&xmesh);

  Z[0]=26;
  Z[1]=34;

  confnum=(int)((gstop-gstart)/gmesh)+1;

  for(iter=0;iter<confnum;iter++) {
    g=gstart+iter*gmesh;
    for(x[1]=xmesh;x[1]<1.0;x[1]+=xmesh) {
      x[0]=1.0-x[1];
      fTCP[0]=f_lTCP(g,x[0],x[1],Z[0],Z[1]);
      fTCP[1]=f_sTCP(g,x[0],x[1],Z[0],Z[1]);
      norm=fabsf(fTCP[1]);
      min=0.0;
      phase=0;
      if(fTCP[0]/norm<min) {
	min=fTCP[0]/norm;
	phase=0;
      } 
      if(fTCP[1]/norm<min) {
	min=fTCP[1]/norm;
	phase=1;
      }
      for(a[1]=xmesh;a[1]<x[1];a[1]+=xmesh) {
	a[0]=1.0-a[1];
	for(A=xmesh;A<1.0;A+=xmesh) {
	  b[1]=(x[1]-A*a[1])/(1.0-A);
	  b[0]=1.0-b[1];
	  fTCP[2]=A*f_lTCP(g,a[0],a[1],Z[0],Z[1])+(1.0-A)*f_sTCP(g,b[0],b[1],Z[0],Z[1]);
	  if(fTCP[2]/norm<min) {
	    min=fTCP[2]/norm;
	    phase=2;
	  }
	}
      }
      //if(phase==0)
      //  printf("C %e %4.2e\n",175/g,x[1]);
      //else if(phase==1)
      //  printf("C %e %4.2e\n",175/g,x[1]);
      //else if(phase==2)
      if(phase==2)
	printf("%e %e\n",175/g,x[1]);
    }
  }
}

double f_l(double g) {
  double value;

  value=-0.899172*g+1.8645*pow(g,0.32301)-0.2748*log(g)-1.4019;

  return value;
}

double f_s(double g) {
  double value;

  value=-0.895929*g+1.5*log(g)-1.1703-10.84/g;

  return value;
}

double Delta_g(double x, double rz) {
  double num[2], denom[2], value;

  num[0]=0.05*(rz-1)*(rz-1);
  denom[0]=(1+0.64*(rz-1))*(1+0.5*(rz-1)*(rz-1));
  num[1]=num[0]/denom[0];
  denom[1]=1+27*(rz-1)*sqrt(x)*(sqrt(x)-0.3)*(sqrt(x)-0.7)*(sqrt(x)-1.0)/(1+0.1*(rz-1));
  value=num[1]/denom[1];

  return value;
}

double f_lTCP(double g, double x1, double x2, double Z1, double Z2) {
  double value, part[2];

  part[0]=f_l(g)+log(x1*Z1/(x1*Z1+x2*Z2));
  part[1]=f_l(g*pow(Z2/Z1,5.0/3.0))+log(x2*Z2/(x1*Z1+x2*Z2));
  value=x1*part[0]+x2*part[1];

  return value;
}

double f_sTCP(double g, double x1, double x2, double Z1, double Z2) {
  double value, part[2];

  part[0]=f_s(g)+log(x1*Z1/(x1*Z1+x2*Z2));
  part[1]=f_s(g*pow(Z2/Z1,5.0/3.0))+log(x2*Z2/(x1*Z1+x2*Z2));
  value=x1*part[0]+x2*part[1]+g*x1*x2*Delta_g(x2,Z2/Z1);

  return value;
}
