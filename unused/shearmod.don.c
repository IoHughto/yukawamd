#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define ALPHA 1.0/137.036
#define HBARC 197.33

void open_xv4a_(char*,int*,float*,float*,float*,int);
void read_xv4a_(char*,int*,float*,float*,float*,int);
void close_xv4a_(char*,int*,float*,float*,float*,int);

int main(int argc, char *argv[])
{
  int j, n, k, PART, a[3], b, c, confnum, conf;
  float time4;
  double BOXSIZE[3], epsilon, PE[6][5], potential, r, relpos[5][3], CHARGE, SCRLEN, CUTOFF, u[3][5][3], v[3][5][3], ahc, first, second, rp[3], rho, time8;
  char line[300], dump[5], filename[100];

  if(argc != 8) {
    fprintf(stderr,"Syntax Error: %s <input file> <# of particles> <density in 1/fm^3> <charge> <epsilon> <screening length> <# of configurations>\n",argv[0]);
    return 0;
  }

  PART=atoi(argv[2]);
  sscanf(argv[3],"%lf",&rho);
  sscanf(argv[4],"%lf",&CHARGE);
  sscanf(argv[5],"%lf",&epsilon);
  sscanf(argv[6],"%lf",&SCRLEN);
  confnum=atoi(argv[7]);
  
  //  double pos[PART][3], vel[PART][3];
  float pos[PART][3], vel[PART][3];

  ahc=ALPHA*HBARC;

  sprintf(filename,"%s",argv[1]);
  k=strlen(filename);
  open_xv4a_(filename,&PART,&time4,&pos[0][0],&vel[0][0],k);
  
  BOXSIZE[0]=pow(PART/rho,1.0/3.0);
  BOXSIZE[1]=BOXSIZE[0];
  BOXSIZE[2]=BOXSIZE[0];

  u[0][0][0]=1.0;
  u[0][0][1]=1.0;
  u[0][0][2]=1.0;
  u[0][1][0]=1+(1+3.0*epsilon/4.0)*epsilon;
  u[0][1][1]=1/sqrt(u[0][1][0]);
  u[0][1][2]=u[0][1][1];
  u[0][2][0]=1-(1-3.0*epsilon/4.0)*epsilon;
  u[0][2][1]=1/sqrt(u[0][2][0]);
  u[0][2][2]=u[0][2][1];
  u[0][3][0]=1+(2.0+3.0*epsilon)*epsilon;
  u[0][3][1]=1/sqrt(u[0][3][0]);
  u[0][3][2]=u[0][3][1];
  u[0][4][0]=1-(2.0-3.0*epsilon)*epsilon;
  u[0][4][1]=1/sqrt(u[0][4][0]);
  u[0][4][2]=u[0][4][1];
  for(b=0;b<5;b++) {
    u[1][b][0]=u[0][b][2];
    u[1][b][1]=u[0][b][0];
    u[1][b][2]=u[0][b][1];
    u[2][b][0]=u[0][b][1];
    u[2][b][1]=u[0][b][2];
    u[2][b][2]=u[0][b][0];
  }
  v[0][0][0]=0.0;
  v[0][0][1]=0.0;
  v[0][0][2]=0.0;
  v[0][1][0]=epsilon/2.0;
  v[0][1][1]=epsilon/2.0;
  v[0][1][2]=epsilon*epsilon/4.0;
  v[0][2][0]=-v[0][1][0];
  v[0][2][1]=-v[0][1][1];
  v[0][2][2]=v[0][1][2];
  v[0][3][0]=2*v[0][1][0];
  v[0][3][1]=2*v[0][1][1];
  v[0][3][2]=4*v[0][1][2];
  v[0][4][0]=-v[0][3][0];
  v[0][4][1]=-v[0][3][1];
  v[0][4][2]=v[0][3][2];
  for(b=0;b<5;b++) {
    v[1][b][0]=v[0][b][2];
    v[1][b][1]=v[0][b][0];
    v[1][b][2]=v[0][b][1];
    v[2][b][0]=v[0][b][1];
    v[2][b][1]=v[0][b][2];
    v[2][b][2]=v[0][b][0];
  }

  for(conf=0;conf<confnum;conf++) {
    sprintf(filename,"%s",argv[1]);
    k=strlen(filename);
    read_xv4a_(filename,&PART,&time4,&pos[0][0],&vel[0][0],k);
    for(b=0;b<6;b++)
      for(c=0;c<5;c++)
	PE[b][c]=0;
    
    for(j=0;j<PART;j++) {
      for(k=0;k<j;k++) {
	if(k!=j) {
	  for(a[0]=0;a[0]<3;a[0]++) {
	    for(a[1]=0;a[1]<3;a[1]++) {
	      for(a[2]=0;a[2]<3;a[2]++) {
		for(n=0;n<3;n++) {
		  rp[n]=pos[j][n]-(pos[k][n]+(a[n]-1)*BOXSIZE[n]);
		  relpos[0][n]=rp[n]*u[0][0][n];
		  relpos[1][n]=rp[n]*u[0][1][n];
		  relpos[2][n]=rp[n]*u[0][2][n];
		  relpos[3][n]=rp[n]*u[0][3][n];
		  relpos[4][n]=rp[n]*u[0][4][n];
		}
		for(n=0;n<5;n++) {
		  r=sqrt(relpos[n][0]*relpos[n][0]+relpos[n][1]*relpos[n][1]+relpos[n][2]*relpos[n][2]);
		  potential=ahc*CHARGE*CHARGE*(exp(-r/SCRLEN)/r);//-exp(-CUTOFF)/(CUTOFF));
		  PE[0][n]+=potential;
		}
		for(n=0;n<3;n++) {
		  relpos[0][n]=rp[n]*u[1][0][n];
		  relpos[1][n]=rp[n]*u[1][1][n];
		  relpos[2][n]=rp[n]*u[1][2][n];
		  relpos[3][n]=rp[n]*u[1][3][n];
		  relpos[4][n]=rp[n]*u[1][4][n];
		}
		for(n=0;n<5;n++) {
		  r=sqrt(relpos[n][0]*relpos[n][0]+relpos[n][1]*relpos[n][1]+relpos[n][2]*relpos[n][2]);
		  potential=ahc*CHARGE*CHARGE*(exp(-r/SCRLEN)/r);//-exp(-CUTOFF)/(CUTOFF));
		  PE[1][n]+=potential;
		}
		for(n=0;n<3;n++) {
		  relpos[0][n]=rp[n]*u[2][0][n];
		  relpos[1][n]=rp[n]*u[2][1][n];
		  relpos[2][n]=rp[n]*u[2][2][n];
		  relpos[3][n]=rp[n]*u[2][3][n];
		  relpos[4][n]=rp[n]*u[2][4][n];
		}
		for(n=0;n<5;n++) {
		  r=sqrt(relpos[n][0]*relpos[n][0]+relpos[n][1]*relpos[n][1]+relpos[n][2]*relpos[n][2]);
		  potential=ahc*CHARGE*CHARGE*(exp(-r/SCRLEN)/r);//-exp(-CUTOFF)/(CUTOFF));
		  PE[2][n]+=potential;
		}
		for(n=1;n<5;n++) {
		  relpos[n][0]=rp[0]+v[0][n][1]*rp[1];
		  relpos[n][1]=rp[1]+v[0][n][0]*rp[0];
		  relpos[n][2]=rp[2]*(1+v[0][n][2]);
		  r=sqrt(relpos[n][0]*relpos[n][0]+relpos[n][1]*relpos[n][1]+relpos[n][2]*relpos[n][2]);
		  potential=ahc*CHARGE*CHARGE*(exp(-r/SCRLEN)/r);//-exp(-CUTOFF)/(CUTOFF));
		  PE[3][n]+=potential;
		}
		for(n=1;n<5;n++) {
		  relpos[n][2]=rp[2]+v[1][n][1]*rp[1];
		  relpos[n][1]=rp[1]+v[1][n][2]*rp[2];
		  relpos[n][0]=rp[0]*(1+v[1][n][0]);
		  r=sqrt(relpos[n][0]*relpos[n][0]+relpos[n][1]*relpos[n][1]+relpos[n][2]*relpos[n][2]);
		  potential=ahc*CHARGE*CHARGE*(exp(-r/SCRLEN)/r);//-exp(-CUTOFF)/(CUTOFF));
		  PE[4][n]+=potential;
		}
		for(n=1;n<5;n++) {
		  relpos[n][2]=rp[2]+v[2][n][0]*rp[0];
		  relpos[n][0]=rp[0]+v[2][n][2]*rp[2];
		  relpos[n][1]=rp[1]*(1+v[2][n][1]);
		  r=sqrt(relpos[n][0]*relpos[n][0]+relpos[n][1]*relpos[n][1]+relpos[n][2]*relpos[n][2]);
		  potential=ahc*CHARGE*CHARGE*(exp(-r/SCRLEN)/r);//-exp(-CUTOFF)/(CUTOFF));
		  PE[5][n]+=potential;
		}
	      }
	    }
	  }
	}
      }
    }
    PE[3][0]=PE[0][0];
    PE[4][0]=PE[1][0];
    PE[5][0]=PE[1][0];
    for(b=0;b<5;b++) {
      PE[0][b]/=PART;
      PE[1][b]/=PART;
      PE[2][b]/=PART;
      PE[3][b]/=PART;
      PE[4][b]/=PART;
      PE[5][b]/=PART;
    }
    printf("%4.3f %16.14e ",epsilon,PE[0][0]);
    for(b=0;b<6;b++) {
      first=(2.0*(PE[b][1]-PE[b][2])/3.0-(PE[b][3]-PE[b][4])/12.0)/epsilon;
      second=(4.0*(PE[b][1]+PE[b][2])/3.0-(PE[b][3]+PE[b][4])/12.0-5.0*PE[b][0]/2.0)/(epsilon*epsilon);
      printf("%16.14e %16.14e ",first,second);
    }
    printf("\n");
  }
  sprintf(filename,"%s",argv[1]);
  k=strlen(filename);
  close_xv4a_(filename,&PART,&time4,&pos[0][0],&vel[0][0],k);
  return 0;
}
