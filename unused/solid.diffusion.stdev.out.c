#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define min(x1,x2) ((x1) > (x2) ? (x2):(x1))

void open_xv4a_(char*,int*,float*,float*,float*,int);
void read_xv4a_(char*,int*,float*,float*,float*,int);
void close_xv4a_(char*,int*,float*,float*,float*,int);
double findshort(double a, double b, double *lr, double BOXSIZE);

void read_xv8_(char*,int*,double*,double*,double*,int);

int main(int argc, char *argv[]) {

  int PART, binnum, k, j, n, flag, coords, confnum, conf, ionnum[5], dumpi;
  float time4, dumpf;
  double time8;
  double maxd, bind, r, x, y, z, diffusion[5], fastd, startt, stopt, cadence, dt, dfac;
  char filename[100], line[300];
  FILE *in, *finalin, *ionfile, *out[5];

  if(argc != 12) {
    fprintf(stderr,"Syntax Error: %s <ion.dat file> <output base> <# of particles> <# of bins> <max distance> <r, x, y, z, or rho> <start time> <stop time> <Delta t> <data cadence> <Delta r threshold>\n",argv[0]);
    return 0;
  }

  ionfile=fopen(argv[1],"r");
  PART=atoi(argv[3]);
  binnum=atoi(argv[4]);
  sscanf(argv[5],"%lf",&maxd);
  coords=atoi(argv[6]);
  sscanf(argv[7],"%lf",&startt);
  sscanf(argv[8],"%lf",&stopt);
  sscanf(argv[9],"%lf",&dt);
  sscanf(argv[10],"%lf",&cadence);
  sscanf(argv[11],"%lf",&fastd);

  switch(coords) {
  case 0:
    dfac=6.0;
    break;
  case 1:
    dfac=2.0;
    break;
  case 2:
    dfac=2.0;
    break;
  case 3:
    dfac=2.0;
    break;
  case 4:
    dfac=4.0;
    break;
  }

  double pos[2][3*PART], skip[PART][3];
  double charge[PART], mass[PART];
  int bin[5][binnum];

  bind=maxd/binnum;
  
  confnum=(int)((stopt-startt-dt)/cadence)+1;
  
  double diff[5][confnum], diffmean[5], stdev[5];

  for(j=0;j<5;j++)
    ionnum[j]=0;

  fgets(line,300,ionfile);

  for(j=0;j<PART;j++) {
    fgets(line,300,ionfile);
    sscanf(line,"%d %lf %lf",&dumpi,&charge[j],&mass[j]);
    if(charge[j]<3.0)
      ionnum[0]++;
    else if(charge[j]<7.0)
      ionnum[1]++;
    else if(charge[j]<9.0)
      ionnum[2]++;
    else if(charge[j]<11.0)
      ionnum[3]++;
    else 
      ionnum[4]++;
  }

  diffusion[0]=0.0; 
  diffusion[1]=0.0; 
  diffusion[2]=0.0; 
  diffusion[3]=0.0; 
  diffusion[4]=0.0;

  for(n=0;n<binnum;n++) {
    bin[0][n]=0;
    bin[1][n]=0;
    bin[2][n]=0;
    bin[3][n]=0;
    bin[4][n]=0;
  }

  for(conf=0;conf<confnum;conf++) {
    for(j=0;j<5;j++)
      diff[j][conf]=0.0;
    //    printf("%d %d\n",conf,confnum);
    sprintf(filename,"md.ckpt.%011.0lf",startt+conf*cadence);
    k=strlen(filename);
    read_xv8_(filename,&PART,&time8,&pos[0][0],&skip[0][0],k);
    sprintf(filename,"md.ckpt.%011.0lf",startt+conf*cadence+dt);
    k=strlen(filename);
    read_xv8_(filename,&PART,&time8,&pos[1][0],&skip[0][0],k);
    
    for(j=0;j<PART;j++) {
      x=min(min(fabsf(pos[1][3*j]-pos[0][3*j]),fabsf(pos[1][3*j]-pos[0][3*j]-maxd)),fabsf(pos[1][3*j]-pos[0][3*j]+maxd));
      y=min(min(fabsf(pos[1][3*j+1]-pos[0][3*j+1]),fabsf(pos[1][3*j+1]-pos[0][3*j+1]-maxd)),fabsf(pos[1][3*j+1]-pos[0][3*j+1]+maxd));
      z=min(min(fabsf(pos[1][3*j+2]-pos[0][3*j+2]),fabsf(pos[1][3*j+2]-pos[0][3*j+2]-maxd)),fabsf(pos[1][3*j+2]-pos[0][3*j+2]+maxd));
      r=sqrt(x*x+y*y+z*z);
      if(coords==1)
	r=x;
      else if(coords==2)
	r=y;
      else if(coords==3)
	r=z;
      else if(coords==4)
	r=sqrt(x*x+y*y);
      flag=0;
      if(r>fastd) {
	if(charge[j]<3.0){
	  diffusion[0]+=r*r;
	  diff[0][conf]+=r*r;
	}
	else if(charge[j]<7.0){
	  diffusion[1]+=r*r;
	  diff[1][conf]+=r*r;
	}
	else if(charge[j]<9.0){
	  diffusion[2]+=r*r;
	  diff[2][conf]+=r*r;
	}
	else if(charge[j]<11.0){
	  diffusion[3]+=r*r;
	  diff[3][conf]+=r*r;
	}
	else{
	  diffusion[4]+=r*r;
	  diff[4][conf]+=r*r;
	}
      }
      for(n=0;n<binnum;n++) {
	if(flag==0) {
	  if(r<(n+1)*bind) {
	    if(charge[j]<3.0)
	      bin[0][n]++;
	    else if(charge[j]<7.0)
	      bin[1][n]++;
	    else if(charge[j]<9.0)
	      bin[2][n]++;
	    else if(charge[j]<11.0)
	      bin[3][n]++;
	    else
	      bin[4][n]++;
	    flag=1;
	  }
	}
      }
    }
  } 
  
  if(ionnum[0]>0) {
    sprintf(filename,"%s.He.out",argv[2]);
    out[0]=fopen(filename,"w");
    for(n=0;n<binnum;n++) {
      fprintf(out[0],"%lf %e\n",n*bind+bind/2.0,(double)bin[0][n]/confnum);
    }
    diffmean[0]=diffusion[0]/(dfac*dt*ionnum[0]*confnum);
    stdev[0]=0.0;
    for(j=0;j<confnum;j++) 
      stdev[0]+=(diffmean[0]-diff[0][j]/(dfac*dt*ionnum[0]))*(diffmean[0]-diff[0][j]/(dfac*dt*ionnum[0]));
    stdev[0]=sqrt(stdev[0]/(confnum-1));
    printf("D_He = %e +/- %e\n",diffmean[0],stdev[0]);
    for(conf=0;conf<confnum;conf++)
      printf("%11.0lf %16.14e\n",startt+conf*cadence,diff[0][conf]/(dfac*dt*ionnum[0]));
  }
  if(ionnum[1]>0) {
    sprintf(filename,"%s.C.out",argv[2]);
    out[1]=fopen(filename,"w");
    for(n=0;n<binnum;n++) {
      fprintf(out[1],"%lf %e\n",n*bind+bind/2.0,(double)bin[1][n]/confnum);
    }
    diffmean[1]=diffusion[1]/(dfac*dt*ionnum[1]*confnum);
    stdev[1]=0.0;
    for(j=0;j<confnum;j++) 
      stdev[1]+=(diffmean[1]-diff[1][j]/(dfac*dt*ionnum[1]))*(diffmean[1]-diff[1][j]/(dfac*dt*ionnum[1]));
    stdev[1]=sqrt(stdev[1]/(confnum-1));
    printf("D_C = %e +/- %e\n",diffmean[1],stdev[1]);
    for(conf=0;conf<confnum;conf++)
      printf("%11.0lf %16.14e\n",startt+conf*cadence,diff[1][conf]/(dfac*dt*ionnum[1]));
  }
  if(ionnum[2]>0) {
    sprintf(filename,"%s.O.out",argv[2]);
    out[2]=fopen(filename,"w");
    for(n=0;n<binnum;n++) {
      fprintf(out[2],"%lf %e\n",n*bind+bind/2.0,(double)bin[2][n]/confnum);
    }
    diffmean[2]=diffusion[2]/(dfac*dt*ionnum[2]*confnum);
    stdev[2]=0.0;
    for(j=0;j<confnum;j++) 
      stdev[2]+=(diffmean[2]-diff[2][j]/(dfac*dt*ionnum[2]))*(diffmean[2]-diff[2][j]/(dfac*dt*ionnum[2]));
    stdev[2]=sqrt(stdev[2]/(confnum-1));
    printf("D_O = %e +/- %e\n",diffmean[2],stdev[2]);
    for(conf=0;conf<confnum;conf++)
      printf("%11.0lf %16.14e\n",startt+conf*cadence,diff[2][conf]/(dfac*dt*ionnum[2]));
  }
  if(ionnum[3]>0) {
    sprintf(filename,"%s.Ne.out",argv[2]);
    out[3]=fopen(filename,"w");
    for(n=0;n<binnum;n++) {
      fprintf(out[3],"%lf %e\n",n*bind+bind/2.0,(double)bin[3][n]/confnum);
    }
    diffmean[3]=diffusion[3]/(dfac*dt*ionnum[3]*confnum);
    stdev[3]=0.0;
    for(j=0;j<confnum;j++) 
      stdev[3]+=(diffmean[3]-diff[3][j]/(dfac*dt*ionnum[3]))*(diffmean[3]-diff[3][j]/(dfac*dt*ionnum[3]));
    stdev[3]=sqrt(stdev[3]/(confnum-1));
    printf("D_Ne = %e +/- %e\n",diffmean[3],stdev[3]);
    for(conf=0;conf<confnum;conf++)
      printf("%11.0lf %16.14e\n",startt+conf*cadence,diff[3][conf]/(dfac*dt*ionnum[3]));
  }
  if(ionnum[4]>0) {
    sprintf(filename,"%s.Fe.out",argv[2]);
    out[4]=fopen(filename,"w");
    for(n=0;n<binnum;n++) {
      fprintf(out[4],"%lf %e\n",n*bind+bind/2.0,(double)bin[4][n]/confnum);
    }
    diffmean[4]=diffusion[4]/(dfac*dt*ionnum[4]*confnum);
    stdev[4]=0.0;
    for(j=0;j<confnum;j++) 
      stdev[4]+=(diffmean[4]-diff[4][j]/(dfac*dt*ionnum[4]))*(diffmean[4]-diff[4][j]/(dfac*dt*ionnum[4]));
    stdev[4]=sqrt(stdev[4]/(confnum-1));
    printf("D_Fe = %e +/- %e\n",diffmean[4],stdev[4]);
    for(conf=0;conf<confnum;conf++)
      printf("%11.0lf %16.14e\n",startt+conf*cadence,diff[4][conf]/(dfac*dt*ionnum[4]));
  }
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
