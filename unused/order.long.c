#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define PI 3.141592653589

void read_xv8_(char*,int*,double*,double*,double*,int);
double findshort(double a, double b, double *lr, double BOXSIZE);
double plgndr(int l, int m, double x);
double ylm(int l, int m, double x);

int main(int argc, char *argv[]) {
  int ii, j, k, l, m, n, mm, PART, flag, dumpi, nn, maxmoment, conf, confnum;
  double theta, phi, rho, BOXSIZE, rp[3], r, time8, lr[3], nncut;
  long long startt, stopt, cadence;
  char line[300], filename[100];
  FILE *in, *ionfile;

  if( argc != 8) {
    fprintf(stderr,"Syntax Error: %s <ion.dat file> <density> <start time> <end time> <data cadence> <nearest neighbor cutoff distance> <max Y_lm moment>\n",argv[0]);
    return 1;
  }

  ionfile=fopen(argv[1],"r");
  sscanf(argv[2],"%lf",&rho);
  startt=atol(argv[3]);
  stopt=atol(argv[4]);
  cadence=atol(argv[5]);
  sscanf(argv[6],"%lf",&nncut);
  maxmoment=atol(argv[7]);

  confnum=(stopt-startt)/cadence+1;

  fgets(line,300,ionfile);
  sscanf(line,"%d",&PART);

  double charge[PART], mass[PART], pos[3*PART], vel[3*PART];
  BOXSIZE=pow(PART/rho,1.0/3.0); 
  double qlm[2][maxmoment+1][2*maxmoment+1], ql[maxmoment+1], value, sum;
  for(j=0;j<PART;j++) {
    fgets(line,300,ionfile);
    sscanf(line,"%d %lf %lf",&dumpi,&charge[j],&mass[j]);
  }
  for(conf=0;conf<confnum;conf++) {
    sprintf(filename,"md.ckpt.%011ld",startt+conf*cadence);
    k=strlen(filename);
    //  printf("%s\n",filename);
    read_xv8_(filename,&PART,&time8,&pos[0],&vel[0],k);

    nn=0;

    for(l=0;l<=maxmoment;l++) {
      ql[l]=0.0;
      for(m=-l;m<=l;m++) {
	qlm[0][l][m+maxmoment]=0.0;
	qlm[1][l][m+maxmoment]=0.0;
      }
    }

    for(j=0;j<PART;j++) {
      for(k=0;k<PART;k++) {
	if(j!=k) {
	  for(n=0;n<3;n++){
	    rp[n]=findshort(pos[3*j+n],pos[3*k+n],&lr[n],BOXSIZE);
	    rp[n]*=lr[n];
	  }
	  r=sqrt(rp[0]*rp[0]+rp[1]*rp[1]+rp[2]*rp[2]);
	  if(r<nncut) {
	    nn++;
	    theta=rp[2]/r;
	    theta=acos(theta);
	    phi=rp[0]/sqrt(rp[0]*rp[0]+rp[1]*rp[1]);
	    phi=acos(phi);
	    if(rp[1]<0) {
	      phi*=-1.0;
	      phi+=2*PI;
	    }
	    for(l=0;l<=maxmoment;l++) {
	      for(m=-l;m<=l;m++) {
		flag=0;
		value=cos(theta);
		if(m<0) {
		  flag=1;
		  mm=-1*m;
		}
		else
		  mm=m;
		value=plgndr(l,mm,value);
		value=ylm(l,mm,value); 
		if(flag==1)
		  value*=pow(-1.0,mm);
		qlm[0][l][m+maxmoment]+=value*cos(m*phi); 
		qlm[1][l][m+maxmoment]+=value*sin(m*phi); 
	      }
	    }
	  }
	}
      }
    }

    printf("%d ",startt+conf*cadence);

    for(l=0;l<=maxmoment;l++) {
      for(m=-l;m<=l;m++) 
	ql[l]+=(qlm[0][l][m+maxmoment]*qlm[0][l][m+maxmoment]+qlm[1][l][m+maxmoment]*qlm[1][l][m+maxmoment])/nn/nn;
      ql[l]*=4*PI/(2*l+1);
      ql[l]=sqrt(ql[l]);
      printf("%e ",l,ql[l]);
    }
    printf("\n");
  }
}

double plgndr(int l, int m, double x) {
  void nrerror(char error_text[]);
  double fact,pll,pmm,pmmp1,somx2;
  int i,ll;


  if (m < 0 || m > l || fabs(x) > 1.0) {
    fprintf(stderr,"Error in P_l calculation\n");
    return 1;
  }
  pmm=1.0;
  if (m > 0) {
    somx2=sqrt((1.0-x)*(1.0+x));
    fact=1.0;
    for (i=1;i<=m;i++) {
      pmm *= -fact*somx2;
      fact += 2.0;
    }
  }
  if (l == m)
    return pmm;
  else {
    pmmp1=x*(2*m+1)*pmm;
    if (l == (m+1))
      return pmmp1;
    else {
      for (ll=m+2;ll<=l;ll++) {
	pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
	pmm=pmmp1;
	pmmp1=pll;
      }
      return pll;
    }
  }
}

double ylm(int l, int m, double x) {
  int ii, flag=0;
  double prefactor, fact, value;
  
  prefactor=(2.0*l+1.0)/(4*PI);
  fact=1.0;
  for(ii=l-m;ii>0;ii--)
    fact*=ii;
  prefactor*=fact;
  fact=1.0;
  for(ii=l+m;ii>0;ii--)
    fact*=ii;
  prefactor/=fact;

  value=sqrt(prefactor)*x;

  return value;
}

double findshort(double a, double b, double *lr, double BOXSIZE) {
  double small, sign;
  
  // Say the closest neighbor is in the left box
  small=fabsf(a-(b-BOXSIZE));
  sign=1.0;
  
  // Check to see if it's in the middle box
  if(fabsf(a-b)<=fabsf(small)) {
    small=fabsf(a-b);
    if(small!=0.0)
      sign=(a-b)/small;
    else
      sign=1.0;
  }
  
  // Finally, see if it's in the right box
  if(fabsf(a-(b+BOXSIZE))<=fabsf(small)) {
    small=fabsf(a-(b+BOXSIZE));
    sign=-1.0;
  }
  
  // Give direction away from nearest neighbor
  *lr=sign;
  
  // Return the smallest value
  return small;
}
