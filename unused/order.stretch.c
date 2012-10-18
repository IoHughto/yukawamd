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
  double theta, phi, rho, xl[3], rp[3], r, time, lr[3], nncut, ev, ek, pp[9], px;
  double startt, stopt, cadence;
  char line[300], filename[100], dumpc[100];
  FILE *in, *ionfile;
  float dump;

  if( argc != 8) {
    fprintf(stderr,"Syntax Error: %s <ion.dat file> <density> <start time> <end time> <data cadence> <nearest neighbor cutoff distance> <max Y_lm moment>\n",argv[0]);
    return 1;
  }

  ionfile=fopen(argv[1],"r");
  sscanf(argv[2],"%lf",&rho);
  sscanf(argv[3],"%lf",&startt);
  sscanf(argv[4],"%lf",&stopt);
  sscanf(argv[5],"%lf",&cadence);
  sscanf(argv[6],"%lf",&nncut);
  maxmoment=atol(argv[7]);

  confnum=(int)((stopt-startt)/cadence)+1;

  fgets(line,300,ionfile);
  sscanf(line,"%d",&PART);

  double charge[PART], mass[PART], pos[3*PART], vel[3*PART];
  double qlm[2][maxmoment+1][2*maxmoment+1], ql[maxmoment+1], value, sum;
  for(j=0;j<PART;j++) {
    fgets(line,300,ionfile);
    sscanf(line,"%d %lf %lf",&dumpi,&charge[j],&mass[j]);
  }
  for(conf=0;conf<confnum;conf++) {
    sprintf(filename,"md.%011.0lf.xv8b",startt+conf*cadence);
    k=strlen(filename);  
    in=fopen(filename,"r");

    fread(&dump,sizeof(float),1,in);
    fread(&dumpc,sizeof(char),6,in);
    fread(&dump,sizeof(float),1,in);
    
    fread(&dump,sizeof(float),1,in);
    fread(&dumpc,sizeof(char),10,in);
    fread(&dumpc,sizeof(char),8,in);
    fread(&dump,sizeof(float),1,in);
    
    fread(&dump,sizeof(float),1,in);
    fread(&dumpc,sizeof(char),8,in);
    fread(&dumpc,sizeof(char),10,in);
    fread(&dumpc,sizeof(char),5,in);
    fread(&dump,sizeof(float),1,in);
    
    fread(&dump,sizeof(float),1,in);
    fread(&dumpc,sizeof(char),20,in);
    fread(&dump,sizeof(float),1,in);
    
    fread(&dump,sizeof(float),1,in);
    fread(&time,sizeof(double),1,in);
    fread(&xl[0],sizeof(double),1,in);
    fread(&xl[1],sizeof(double),1,in);
    fread(&xl[2],sizeof(double),1,in);
    fread(&ev,sizeof(double),1,in);
    fread(&ek,sizeof(double),1,in);
    fread(&px,sizeof(double),1,in);
    fread(&pp[0],sizeof(double),9,in);
    fread(&dumpi,sizeof(int),1,in);
    fread(&dump,sizeof(float),1,in);
    
    fread(&dump,sizeof(float),1,in);

    for(j=0;j<PART;j++) {
      fread(&pos[3*j],sizeof(double),1,in);
      fread(&pos[3*j+1],sizeof(double),1,in);
      fread(&pos[3*j+2],sizeof(double),1,in);
      fread(&vel[3*j],sizeof(double),1,in);
      fread(&vel[3*j+1],sizeof(double),1,in);
      fread(&vel[3*j+2],sizeof(double),1,in);
    }

    fread(&dump,sizeof(float),1,in);

    fclose(in);

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
	    rp[n]=findshort(pos[3*j+n],pos[3*k+n],&lr[n],xl[n]);
	    rp[n]*=lr[n];
	  }
	  r=sqrt(rp[0]*rp[0]+rp[1]*rp[1]+rp[2]*rp[2]);
	  if(r<nncut) {
	    nn++;
	    theta=rp[2]/r;
	    theta=acos(theta);
	    if(rp[0]*rp[0]+rp[1]*rp[1]<0.00001)
	      phi=0.0;
	    else
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
