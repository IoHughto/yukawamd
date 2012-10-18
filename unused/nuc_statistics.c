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
  int ii, j, k, l, m, n, mm, PART, flag, dumpi, nn, maxmoment, conf, confnum, comp[2][3], rect, xyzflag;
  double theta, phi, rho, BOXSIZE[3], rp[3], r, time8, lr[3], nncut;
  double startt, stopt, cadence;
  char line[300], filename[100];
  FILE *in, *ionfile, *xyz;

  if( argc != 7) {
    fprintf(stderr,"Syntax Error: %s <ion.dat file> <density> <start time> <end time> <data cadence> <nearest neighbor cutoff distance>\n",argv[0]);
    return 1;
  }

  ionfile=fopen(argv[1],"r");
  sscanf(argv[2],"%lf",&rho);
  sscanf(argv[3],"%lf",&startt);
  sscanf(argv[4],"%lf",&stopt);
  sscanf(argv[5],"%lf",&cadence);
  sscanf(argv[6],"%lf",&nncut);

  confnum=(int)((stopt-startt)/cadence)+1;

  xyzflag=0;

  if(confnum==1) {
    xyzflag=1;
    sprintf(filename,"md.%011.0lf.xyz",startt);
    xyz=fopen(filename,"w");
  }

  fgets(line,300,ionfile);
  sscanf(line,"%d",&PART);

  double charge[PART], mass[PART], pos[3*PART], vel[3*PART];

  if(rect==1) {
    BOXSIZE[0]=pow(PART/(2.0*rho),1.0/3.0); 
    BOXSIZE[1]=pow(PART/(2.0*rho),1.0/3.0); 
    BOXSIZE[2]=2.0*pow(PART/(2.0*rho),1.0/3.0); 
  }
  else {
    BOXSIZE[0]=pow(PART/rho,1.0/3.0);
    BOXSIZE[1]=pow(PART/rho,1.0/3.0);
    BOXSIZE[2]=pow(PART/rho,1.0/3.0);
  }
  int near[PART][28], nearc[PART], bond, nearcomp[PART][3], solid[confnum][PART], interface[PART][500], interfacec[PART], fsolid[PART];
  double qlm[PART][4][13], qlbar[2][PART], qlmbar[PART][4][13], value, sum, qldot;
  for(j=0;j<PART;j++) {
    fgets(line,300,ionfile);
    sscanf(line,"%d %lf %lf",&dumpi,&charge[j],&mass[j]);
  }
  
  for(j=0;j<PART;j++) {
    for(k=0;k<500;k++)
      interface[j][k]=0;
    interfacec[j]=0;
  }

  comp[0][0]=0;
  comp[0][1]=0;
  comp[0][2]=0;
  comp[1][0]=0;
  comp[1][1]=0;
  comp[1][2]=0;
    
  for(conf=0;conf<confnum;conf++) {
 
    for(j=0;j<PART;j++) {
      for(k=0;k<28;k++)
	near[j][k]=0;
      for(k=0;k<3;k++)
	nearcomp[j][k]=0;
      nearc[j]=0;
    }

    sprintf(filename,"md.ckpt.%011.0lf",startt+conf*cadence);
    k=strlen(filename);
    read_xv8_(filename,&PART,&time8,&pos[0],&vel[0],k);

    nn=0;

    for(j=0;j<PART;j++) {   
      nn=0;
      for(m=-6;m<=6;m++) {
	qlm[j][0][m+6]=0.0;
	qlm[j][1][m+6]=0.0;
	qlm[j][2][m+6]=0.0;
	qlm[j][3][m+6]=0.0;
      }
      
      for(k=0;k<PART;k++) {
	if(j!=k) {
	  for(n=0;n<3;n++){
	    rp[n]=findshort(pos[3*j+n],pos[3*k+n],&lr[n],BOXSIZE[n]);
	    rp[n]*=lr[n];
	  }
	  r=sqrt(rp[0]*rp[0]+rp[1]*rp[1]+rp[2]*rp[2]);
	  if(r<50.0) {
	    interface[j][interfacec[j]]=k;
	    interfacec[j]++;
	  }
	  if(r<nncut) {
	    near[j][nearc[j]]=k;
	    if(charge[k]<7.0)
	      nearcomp[j][0]++;
	    else if(charge[k]<9.0)
	      nearcomp[j][1]++;
	    else
	      nearcomp[j][2]++;
	    nearc[j]++;
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
	    for(m=-6;m<=6;m++) {
	      flag=0;
	      value=cos(theta);
	      if(m<0) {
		flag=1;
		mm=-1*m;
	      }
	      else
		mm=m;
	      value=plgndr(6,mm,value);
	      value=ylm(6,mm,value); 
	      if(flag==1)
		value*=pow(-1.0,mm);
	      qlm[j][0][m+6]+=value*cos(m*phi); 
	      qlm[j][1][m+6]+=value*sin(m*phi); 
	    }
	    for(m=-4;m<=4;m++) {
	      flag=0;
	      value=cos(theta);
	      if(m<0) {
		flag=1;
		mm=-1*m;
	      }
	      else
		mm=m;
	      value=plgndr(4,mm,value);
	      value=ylm(4,mm,value); 
	      if(flag==1)
		value*=pow(-1.0,mm);
	      qlm[j][2][m+4]+=value*cos(m*phi); 
	      qlm[j][3][m+4]+=value*sin(m*phi); 
	    }
	  }
	}
      }
      for(m=-6;m<=6;m++) { 
      	qlm[j][0][m+6]/=nearc[j]; 
       	qlm[j][1][m+6]/=nearc[j]; 
      } 
      for(m=-4;m<=4;m++) { 
       	qlm[j][2][m+4]/=nearc[j]; 
       	qlm[j][3][m+4]/=nearc[j]; 
      } 
    }
    for(j=0;j<PART;j++) {
      for(l=0;l<4;l++)
	for(m=0;m<13;m++)
	  qlmbar[j][l][m]=0.0;
      for(k=0;k<nearc[j];k++) {
	for(m=-6;m<=6;m++) {
	  qlmbar[j][0][m+6]+=qlm[k][0][m+6];
	  qlmbar[j][1][m+6]+=qlm[k][1][m+6];
	}
	for(m=-4;m<=4;m++) {
	  qlmbar[j][2][m+4]+=qlm[k][2][m+4];
	  qlmbar[j][3][m+4]+=qlm[k][3][m+4];
	}
      }
      for(m=-6;m<=6;m++) {
	qlmbar[j][0][m+6]/=(1.0+nearc[j]);
	qlmbar[j][1][m+6]/=(1.0+nearc[j]);
      }
      for(m=-4;m<=4;m++) {
	qlmbar[j][2][m+4]/=(1.0+nearc[j]);
	qlmbar[j][3][m+4]/=(1.0+nearc[j]);
      }
    }
    for(j=0;j<PART;j++) {
      qlbar[0][j]=0.0;
      qlbar[1][j]=0.0;
      for(m=-6;m<=6;m++) 
	qlbar[0][j]+=qlmbar[j][0][m+6]*qlmbar[j][0][m+6]+qlmbar[j][1][m+6]*qlmbar[j][1][m+6];
      for(m=-4;m<=4;m++) 
	qlbar[1][j]+=qlmbar[j][2][m+4]*qlmbar[j][2][m+4]+qlmbar[j][3][m+4]*qlmbar[j][3][m+4];
      qlbar[0][j]=sqrt(4.0*3.14159/13.0*qlbar[0][j]);
      qlbar[1][j]=sqrt(4.0*3.14159/9.0*qlbar[1][j]);
      printf("%d %lf %lf\n",j,qlbar[1][j],qlbar[0][j]);
    }
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