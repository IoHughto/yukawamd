#include <stdio.h>
#include <math.h>

#define PI 3.141592653589

double plgndr(int l, int m, double x);
double ylm(int l, int m, double x);

int main(int argc, char *argv[]) {
  int ii, j, l, m, mm, PART, flag;
  double theta, phi;
  char line[300];
  FILE *in;

  if( argc != 3) {
    fprintf(stderr,"Syntax Error: %s <theta/phi file> <# of particles>\n",argv[0]);
    return 1;
  }

  in=fopen(argv[1],"r");
  PART=atoi(argv[2]);

  double qlm[2][15][29], ql[15], value, sum;

  for(l=0;l<15;l++) {
    ql[l]=0.0;
    for(m=-l;m<=l;m++) {
      qlm[0][l][m+14]=0.0;
      qlm[1][l][m+14]=0.0;
    }
  }
  
  for(j=0;j<PART;j++) {
    fgets(line,300,in);
    sscanf(line,"%lf %lf",&theta,&phi);
    for(l=0;l<15;l++)
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
	qlm[0][l][m+14]+=value*cos(m*phi); 
	qlm[1][l][m+14]+=value*sin(m*phi); 
      }
  }
  /* for(l=0;l<15;l++)  */
  /*   for(m=-l;m<=l;m++) */
  /*     printf("%d %d %lf\n",l,m,qlm[l][m+14]); */
  for(l=0;l<15;l++) {
    sum=0.0;
    for(m=-l;m<=l;m++) {
      sum+=qlm[0][l][m+14]/PART;
      ql[l]+=(qlm[0][l][m+14]*qlm[0][l][m+14]+qlm[1][l][m+14]*qlm[1][l][m+14])/PART/PART;
    }
    //    ql[l]=(4*PI)/(2.0*l+1.0)*sum*sum;
    ql[l]*=4*PI/(2*l+1);
    ql[l]=sqrt(ql[l]);
    printf("%d %e\n",l,ql[l]);
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
