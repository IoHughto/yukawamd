#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define min(x1,x2) ((x1) > (x2) ? (x2):(x1))

void open_xv4a_(char*,int*,float*,float*,float*,int);
void read_xv4a_(char*,int*,float*,float*,float*,int);
void close_xv4a_(char*,int*,float*,float*,float*,int);
double findshort(double a, double b, double *lr, double BOXSIZE);
char* getspec(double charge, char *specname);

void read_xv8_(char*,int*,double*,double*,double*,int);

int main(int argc, char *argv[]) {

  int PART, binnum, k, j, n, flag, coords, confnum, conf, dumpi;
  float time4, dumpf;
  double time8;
  double maxd, bind, r, x, y, z, fastd, startt, stopt, cadence, dt, dfac;
  char filename[100], line[300], iname[5];
  FILE *in, *finalin, *ionfile, *out, *tout;

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
  double charge[PART], mass[PART], ionvec[100];
  int ioncount;

  bind=maxd/binnum;
  
  confnum=(int)((stopt-startt-dt)/cadence)+1;

  fgets(line,300,ionfile);

  ioncount=0;

  for(j=0;j<PART;j++) {
    fgets(line,300,ionfile);
    sscanf(line,"%d %lf %lf",&dumpi,&charge[j],&mass[j]);
    if(j==0){
      ionvec[ioncount]=charge[j];
      ioncount++;
    }
    for(k=0;k<ioncount;k++){
      if(charge[j]==ionvec[k])
	break;
      if(k==ioncount-1) {
	ionvec[ioncount]=charge[j];
	ioncount++;
      }
    }
  }

  int ionnum[ioncount], bin[ioncount][binnum], ispec;
  double diffusion[ioncount];
  double diff[ioncount][confnum], diffmean[ioncount], stdev[ioncount];

  for(j=0;j<ioncount;j++) 
    ionnum[j]=0;

  for(j=0;j<PART;j++) 
    for(k=0;k<ioncount;k++) 
      if(charge[j]==ionvec[k]) 
	ionnum[k]++;

  for(j=0;j<ioncount;j++)
    diffusion[j]=0.0;

  for(n=0;n<binnum;n++) 
    for(j=0;j<ioncount;j++)
      bin[j][n]=0;

  for(conf=0;conf<confnum;conf++) {
    for(j=0;j<ioncount;j++)
      diff[j][conf]=0.0;
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
      if(r>fastd) {
	for(k=0;k<ioncount;k++) 
	  if(ionvec[k]==charge[j]) {
	    ispec=k;
	    break;
	  }
	diffusion[ispec]+=r*r;
	diff[ispec][conf]+=r*r;
      }
      for(n=0;n<binnum;n++) {
	if(r<(n+1)*bind) {
	  bin[ispec][n]++;
	  break;
	}
      }
    }
  } 
  
  for(j=0;j<ioncount;j++) {
    sprintf(filename,"%s.%d.out",argv[2],(int)ionvec[j]);
    out=fopen(filename,"w");
    sprintf(filename,"%s.%d.time.out",argv[2],(int)ionvec[j]);
    tout=fopen(filename,"w");
    for(n=0;n<binnum;n++) {
      fprintf(out,"%lf %e\n",n*bind+bind/2.0,(double)bin[j][n]/confnum);
    }
    diffmean[j]=diffusion[j]/(dfac*dt*ionnum[j]*confnum);
    stdev[j]=0.0;
    for(k=0;k<confnum;k++) 
      stdev[j]+=(diffmean[j]-diff[j][k]/(dfac*dt*ionnum[j]))*(diffmean[j]-diff[j][k]/(dfac*dt*ionnum[j]));
    stdev[j]=sqrt(stdev[j]/(confnum-1));
    for(conf=0;conf<confnum;conf++)
      fprintf(tout,"%11.0lf %16.14e\n",startt+conf*cadence,diff[j][conf]/(dfac*dt*ionnum[j]));
    getspec(ionvec[j],iname);
    printf("D_%s = %e +/- %e\n",iname,diffmean[j],stdev[j]);
  }

  return 0;
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

char* getspec(double charge, char *specname) {
  char buf[2];

  if (charge<1.5)
    strcpy(buf,"H");
  else if(charge<2.5)
    strcpy(buf,"He");
  else if(charge<3.5)
    strcpy(buf,"Li");
  else if(charge<4.5)
    strcpy(buf,"Be");
  else if(charge<5.5)
    strcpy(buf,"B");
  else if(charge<6.5)
    strcpy(buf,"C");
  else if(charge<7.5)
    strcpy(buf,"N");
  else if(charge<8.5)
    strcpy(buf,"O");
  else if(charge<9.5)
    strcpy(buf,"F");
  else if(charge<10.5)
    strcpy(buf,"Ne");
  else if(charge<11.5)
    strcpy(buf,"Na");
  else if(charge<12.5)
    strcpy(buf,"Mg");
  else if(charge<13.5)
    strcpy(buf,"Al");
  else if(charge<14.5)
    strcpy(buf,"Si");
  else if(charge<15.5)
    strcpy(buf,"P");
  else if(charge<16.5)
    strcpy(buf,"S");
  else if(charge<17.5)
    strcpy(buf,"Cl");
  else if(charge<18.5)
    strcpy(buf,"Ar");
  else if(charge<19.5)
    strcpy(buf,"K");
  else if(charge<20.5)
    strcpy(buf,"Ca");
  else if(charge<21.5)
    strcpy(buf,"Sc");
  else if(charge<22.5)
    strcpy(buf,"Ti");
  else if(charge<23.5)
    strcpy(buf,"V");
  else if(charge<24.5)
    strcpy(buf,"Cr");
  else if(charge<25.5)
    strcpy(buf,"Mn");
  else if(charge<26.5)
    strcpy(buf,"Fe");
  else if(charge<27.5)
    strcpy(buf,"Co");
  else if(charge<28.5)
    strcpy(buf,"Ni");
  else if(charge<29.5)
    strcpy(buf,"Cu");
  else if(charge<30.5)
    strcpy(buf,"Zn");
  else if(charge<31.5)
    strcpy(buf,"Ga");
  else if(charge<32.5)
    strcpy(buf,"Ge");
  else if(charge<33.5)
    strcpy(buf,"As");
  else if(charge<34.5)
    strcpy(buf,"Se");
  else if(charge<35.5)
    strcpy(buf,"Br");
  else if(charge<36.5)
    strcpy(buf,"Kr");
  else if(charge<37.5)
    strcpy(buf,"Rb");
  else if(charge<38.5)
    strcpy(buf,"Sr");
  else if(charge<39.5)
    strcpy(buf,"Y");
  else if(charge<40.5)
    strcpy(buf,"Zr");
  else if(charge<41.5)
    strcpy(buf,"Nb");
  else if(charge<42.5)
    strcpy(buf,"Mo");
  else if(charge<43.5)
    strcpy(buf,"Tc");
  else if(charge<44.5)
    strcpy(buf,"Ru");
  else if(charge<45.5)
    strcpy(buf,"Rh");
  else if(charge<46.5)
    strcpy(buf,"Pd");
  else if(charge<47.5)
    strcpy(buf,"Ag");
  else if(charge<48.5)
    strcpy(buf,"Cd");
  else if(charge<49.5)
    strcpy(buf,"In");
  else
    strcpy(buf,"Sn");

  return strcpy(specname,buf);
}
