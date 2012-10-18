#include <stdio.h>
#include <math.h>

#define PI 3.141592653589

int main(int argc, char *argv[]) {

  int PART, j, binnum, n, tflag, pflag, peaknum, peak, thresh;
  double bind;
  char line[300], dumpc, filename[100];
  FILE *in, *out, *histout;

  if(argc!=6) {
    fprintf(stderr,"Syntax Error: %s <input file> <output file> <# of bins> <# of peaks> <threshold>\n",argv[0]);
    return 1;
  }

  in=fopen(argv[1],"r");
  out=fopen(argv[2],"w");
  histout=fopen("histo.tp","w");
  binnum=atoi(argv[3]);
  peaknum=atoi(argv[4]);
  thresh=atoi(argv[5]);

  bind=PI/binnum;

  int bins[binnum][2*binnum];

  for(n=0;n<binnum;n++) {
    for(j=0;j<binnum;j++) {
      bins[n][j]=0;
      bins[n][j+binnum]=0;
    }
  }

  fgets(line,300,in);
  sscanf(line,"%d",&PART);
  fgets(line,300,in);
  double pos[3][PART], theta[PART], phi[PART];
  
  for(j=0;j<PART;j++) {
    fgets(line,300,in);
    sscanf(line,"%c %lf %lf %lf",&dumpc,&pos[0][j],&pos[1][j],&pos[2][j]);
    theta[j]=pos[2][j]/sqrt(pos[0][j]*pos[0][j]+pos[1][j]*pos[1][j]+pos[2][j]*pos[2][j]);
    theta[j]=acos(theta[j]);
    phi[j]=pos[0][j]/sqrt(pos[0][j]*pos[0][j]+pos[1][j]*pos[1][j]);
    phi[j]=acos(phi[j]);
    if(pos[1][j]<0) {
      phi[j]*=-1.0;
      phi[j]+=2*PI;
    }
    fprintf(out,"%e %e\n",theta[j],phi[j]);
    tflag=binnum;
    for(n=0;n<binnum;n++) {
      if(tflag==binnum) {
	if(theta[j]<bind*(n+1)){
	  tflag=n;
	}
      } 
    }
    pflag=2*binnum;
    for(n=0;n<2*binnum;n++) {
      if(pflag==2*binnum) {
	if(phi[j]<bind*(n+1)) {
	  pflag=1;
	  bins[tflag][n]++;
	}
      }
    }
  }

  FILE *binfile;
  sprintf(filename,"%s.bins",argv[2]);
  binfile=fopen(filename,"w");
  for(n=0;n<binnum;n++)
    for(j=0;j<2*binnum;j++)
      fprintf(binfile,"%d %d %d\n",n,j,bins[n][j]);

  int max, tmax, pmax, rad, flag, rcut;
  double r, tmin, pmin;
  int ansum;
  double tsum=0.0, psum=0.0;
  int cluster=0;
    
  for(peak=0;peak<peaknum;peak++) {
    max=0;
    // Find first peak
    for(n=0;n<binnum;n++)
      for(j=0;j<2*binnum;j++)
	if(bins[n][j]>max) {
	  max=bins[n][j];
	  tmax=n;
	  pmax=j;
	}
    
    // See where annulus goes to zero
    flag=0;
    for(rad=0;rad<binnum;rad++) {
      if(flag==0){
	ansum=0.0;
	for(n=0;n<binnum;n++) {
	  for(j=0;j<2*binnum;j++) {
	    tmin=tmax-n;
	    pmin=pmax-j;
	    if(abs(pmax-j+2*binnum)<abs(pmin))
	      pmin=pmax-j+2*binnum;
	    if(abs(pmax-j-2*binnum)<abs(pmin))
	      pmin=pmax-j-2*binnum;
	    r=sqrt(tmin*tmin+pmin*pmin);
	    if(r<(rad+2) && r>=rad) 
	      ansum+=bins[n][j];
	  }
	}
	if(ansum<thresh) {
	  flag=1;
	  rcut=rad;
	}
      }
    }
    
    // Average points inside rings and zero out values
    tsum=0.0;
    psum=0.0;
    cluster=0;
    
    for(j=0;j<PART;j++) {
      flag=0;
      tmin=theta[j]-bind*tmax;
      pmin=phi[j]-bind*pmax;
      if(fabsf(phi[j]-bind*pmax+2*PI)<fabsf(pmin)) {
	pmin=phi[j]-bind*pmax+2*PI;
	flag=1;
      }
      if(fabsf(phi[j]-bind*pmax-2*PI)<fabsf(pmin)) {
      	pmin=phi[j]-bind*pmax-2*PI;
	flag=-1;
      }
      r=sqrt(tmin*tmin+pmin*pmin);
      if(r<rcut*bind) {
	tsum+=theta[j];
	psum+=phi[j]+2*PI*flag;
	cluster++;
      }
    }
    tsum/=cluster;
    psum/=cluster;
    if(psum<0)
      psum+=2*PI;
    printf("%e %e\n",tsum,psum);

    // Remove counts from peak
    for(n=0;n<binnum;n++)
      for(j=0;j<2*binnum;j++){
	tmin=tmax-n;
	pmin=pmax-j;
	if(abs(pmax-j+2*binnum)<abs(pmin))
	  pmin=pmax-j+2*binnum;
	if(abs(pmax-j-2*binnum)<abs(pmin))
	  pmin=pmax-j-2*binnum;
	r=sqrt(tmin*tmin+pmin*pmin);
	if(r<rcut)
	  bins[n][j]=0;
      }   
  }
   for(n=0;n<binnum;n++) 
     for(j=0;j<2*binnum;j++) 
       fprintf(histout,"%e %e %d\n",n*bind,j*bind,bins[n][j]); 
}
