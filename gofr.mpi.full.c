#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpi.h>

#define min(x1,x2) ((x1) > (x2) ? (x2):(x1))

void read_xv8_(char*,int*,double*,double*,double*,int);
double findshort(double a, double b, double *lr, double BOXSIZE);
char* getspec(double charge, char *specname);

#define PI 3.141592653589
#define ALPHA 1.0/137.035999
#define HBARC 197.32693

int main(int argc, char *argv[]) {

  MPI_Status status;
  int PART, binnum, n, j, k, l, iter, flag, dt, dumpi, ierr, rank, count, compnum;
  long long startt, stopt, confnum;
  double rho, BOXSIZE, dump, relpos[3], r, bind, time8;
  char line[300], filename[100], jname[3], kname[3];
  FILE *out, *ionfile;

  if(argc != 9) {
    fprintf(stderr,"Syntax Error: %s <ion.dat file> <output base> <# of particles> <density> <# of bins> <start time> <end time> <data cadence>\n",argv[0]);
    return 0;
  }

  ierr=MPI_Init(&argc,&argv);
  ierr=MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  ierr=MPI_Comm_size(MPI_COMM_WORLD,&count);

  ionfile=fopen(argv[1],"r");
  PART=atoi(argv[3]);
  sscanf(argv[4],"%lf",&rho);
  binnum=atoi(argv[5]);
  startt=atol(argv[6]);
  stopt=atol(argv[7]);
  dt=atoi(argv[8]);
  compnum=PART/count;

  confnum=(stopt-startt)/dt;

  BOXSIZE=pow(PART/rho,1.0/3.0);
  bind=BOXSIZE/(binnum*2.0);

  fgets(line,300,ionfile);

  double charge[PART], mass[PART], ionvec[10];
  int ioncount;

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

  int ionnum[ioncount];

  for(j=0;j<ioncount;j++) 
    ionnum[j]=0;

  for(j=0;j<PART;j++) 
    for(k=0;k<ioncount;k++) 
      if(charge[j]==ionvec[k]) 
	ionnum[k]++;

  double vel[3*PART], pos[3*PART]; 

  int bins[ioncount][ioncount][binnum], recv[ioncount][ioncount][binnum],jspec,kspec; 

  for(l=0;l<binnum;l++)  
    for(j=0;j<ioncount;j++) 
      for(k=0;k<ioncount;k++) 
	bins[j][k][l]=0; 
  
  for(iter=0;iter<confnum;iter++) { 
    
    if(rank==0) { 
      sprintf(filename,"md.ckpt.%011ld",startt+dt*iter); 
      printf("%s\n",filename); 
      k=strlen(filename); 
      read_xv8_(filename,&PART,&time8,&pos[0],&vel[0],k); 
    } 

    ierr=MPI_Bcast(pos,3*PART,MPI_DOUBLE,0,MPI_COMM_WORLD); 

    for(j=rank*compnum;j<(rank+1)*compnum;j++) {  
      for(k=0;k<j;k++) { 
	for(l=0;l<ioncount;l++) 
	  if(ionvec[l]==charge[j]) {
	    jspec=l;
	    break;
	  }
	for(l=0;l<ioncount;l++) 
	  if(ionvec[l]==charge[k]) {
	    kspec=l;
	    break;
	  }
	for(n=0;n<3;n++) { 
   	  relpos[n]=findshort(pos[3*j+n],pos[3*k+n],&dump,BOXSIZE); 
   	} 
   	r=sqrt(relpos[0]*relpos[0]+relpos[1]*relpos[1]+relpos[2]*relpos[2]); 
   	flag=0; 
   	for(n=0;n<binnum;n++) { 
	  if(r<(n+1)*bind) { 
	    bins[jspec][kspec][n]++;
	    break; 
	  } 
	} 
      } 
    } 
  } 
  
  ierr=MPI_Allreduce(bins,recv,ioncount*ioncount*binnum,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

  if(rank==0) { 
    for(j=0;j<ioncount;j++) {
      getspec(ionvec[j],jname);
      for(k=0;k<=j;k++) {
	getspec(ionvec[k],kname);
	sprintf(filename,"%s.%s-%s.out",argv[2],jname,kname);
	out=fopen(filename,"w");
	for(n=0;n<binnum;n++) {
	  r=(double)n/(double)binnum*BOXSIZE/2.0;
	  if(j==k)
	    fprintf(out,"%lf %e\n",r,(double)recv[j][k][n]*PART/(2*PI*r*r*rho*bind*ionnum[j]*ionnum[k]*confnum));
	  else
	    fprintf(out,"%lf %e\n",r,((double)recv[j][k][n]+recv[k][j][n])*PART/(4*PI*r*r*rho*bind*ionnum[j]*ionnum[k]*confnum));
	}
	fclose(out);
      }
    }
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

