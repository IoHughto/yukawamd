#include <stdio.h>
#include <string.h>
#include <math.h>

#define AHC 197.33/137.036

void open_xv4a_(char*,int*,float*,float*,float*,int);
void read_xv4a_(char*,int*,float*,float*,float*,int);
void close_xv4a_(char*,int*,float*,float*,float*,int);
double comppot(int class, double dx, int axis, float* pos, double* charge);

int PART;
double lambda, BOXSIZE;

int main(int argc, char *argv[]) {
  
  int j, k, kk, n, PART, confnum, dumpi, iter[2], group;
  double rho, epsilon, T, vol, pot[6][5], d[2][6], f[3][6], g[3][6],h[5][6], mu[5];
  char line[300], filename[100];
  FILE *dat;

  if (argc != 9) {
    fprintf(stderr,"Syntax Error: %s <input file> <charge file> <# of particles> <# of configs> <density(fm^-3)> <epsilon> <screening length(fm)> <T(MeV)>\n",argv[0]);
  }

  dat=fopen(argv[2],"r");
  PART=atoi(argv[3]);
  confnum=atoi(argv[4]);
  sscanf(argv[5],"%lf",&rho);
  sscanf(argv[6],"%lf",&epsilon);
  sscanf(argv[7],"%lf",&lambda);
  sscanf(argv[8],"%lf",&T);

  group=confnum/4;
  vol=(double)PART/rho;
  BOXSIZE=pow(vol,1.0/3.0);

  double charge[PART], mass[PART];
  float pos[3*PART], vel[3*PART], time4;
 
  fgets(line,300,dat);

  for(j=0;j<PART;j++) {
    fgets(line,300,dat);
    sscanf(line,"%d %lf %lf\n",&dumpi,&charge[j],&mass[j]);
  }
    
  for(j=0;j<6;j++)
    for(k=0;k<3;k++)
      g[k][j]=0.0;
    
  for(j=0;j<6;j++)
    for(k=0;k<5;k++)
      h[k][j]=0.0;
  
  sprintf(filename,"%s",argv[1]);
  kk=strlen(filename);
  open_xv4a_(filename,&PART,&time4,&pos[0],&vel[0],kk);

  for(iter[0]=0;iter[0]<4;iter[0]++) {  
    for(j=0;j<6;j++)
      for(k=0;k<3;k++)
	f[k][j]=0.0;
    for(iter[1]=0;iter[1]<group;iter[1]++) {
      sprintf(filename,"%s",argv[1]);
      kk=strlen(filename);
      read_xv4a_(filename,&PART,&time4,&pos[0],&vel[0],kk);
      
      // Distortion one
      pot[0][0]=comppot(0,-2.0*epsilon,0,&pos[0],&charge[0]);
      pot[0][1]=comppot(0,-1.0*epsilon,0,&pos[0],&charge[0]);
      pot[0][2]=comppot(0,0.0,0,&pos[0],&charge[0]);
      pot[0][3]=comppot(0,epsilon,0,&pos[0],&charge[0]);
      pot[0][4]=comppot(0,2.0*epsilon,0,&pos[0],&charge[0]);

      // Distortion two
      pot[1][0]=comppot(0,-2.0*epsilon,1,&pos[0],&charge[0]);
      pot[1][1]=comppot(0,-1.0*epsilon,1,&pos[0],&charge[0]);
      pot[1][2]=comppot(0,0.0,1,&pos[0],&charge[0]);
      pot[1][3]=comppot(0,epsilon,1,&pos[0],&charge[0]);
      pot[1][4]=comppot(0,2.0*epsilon,1,&pos[0],&charge[0]);

      // Distortion three
      pot[2][0]=comppot(0,-2.0*epsilon,2,&pos[0],&charge[0]);
      pot[2][1]=comppot(0,-1.0*epsilon,2,&pos[0],&charge[0]);
      pot[2][2]=comppot(0,0.0,2,&pos[0],&charge[0]);
      pot[2][3]=comppot(0,epsilon,2,&pos[0],&charge[0]);
      pot[2][4]=comppot(0,2.0*epsilon,2,&pos[0],&charge[0]);

      // Distortion four
      pot[3][0]=comppot(1,-2.0*epsilon,0,&pos[0],&charge[0]);
      pot[3][1]=comppot(1,-1.0*epsilon,0,&pos[0],&charge[0]);
      pot[3][2]=comppot(1,0.0,0,&pos[0],&charge[0]);
      pot[3][3]=comppot(1,epsilon,0,&pos[0],&charge[0]);
      pot[3][4]=comppot(1,2.0*epsilon,0,&pos[0],&charge[0]);

      // Distortion five
      pot[4][0]=comppot(1,-2.0*epsilon,1,&pos[0],&charge[0]);
      pot[4][1]=comppot(1,-1.0*epsilon,1,&pos[0],&charge[0]);
      pot[4][2]=comppot(1,0.0,1,&pos[0],&charge[0]);
      pot[4][3]=comppot(1,epsilon,1,&pos[0],&charge[0]);
      pot[4][4]=comppot(1,2.0*epsilon,1,&pos[0],&charge[0]);

      // Distortion six
      pot[5][0]=comppot(1,-2.0*epsilon,2,&pos[0],&charge[0]);
      pot[5][1]=comppot(1,-1.0*epsilon,2,&pos[0],&charge[0]);
      pot[5][2]=comppot(1,0.0,2,&pos[0],&charge[0]);
      pot[5][3]=comppot(1,epsilon,2,&pos[0],&charge[0]);
      pot[5][4]=comppot(1,2.0*epsilon,2,&pos[0],&charge[0]);

      for(j=0;j<6;j++) {
	d[0][j]=(4.0/3.0*(pot[j][3]+pot[j][1])-(pot[j][4]+pot[j][0])/12.0-2.5*pot[j][2])/(epsilon*epsilon);
	d[1][j]=(2.0/3.0*(pot[j][3]-pot[j][1])-(pot[j][4]-pot[j][0])/12.0)/(epsilon);
	f[0][j]+=d[0][j];
	f[1][j]+=d[1][j];
	f[2][j]+=d[1][j]*d[1][j];
	g[0][j]+=d[0][j];
	g[1][j]+=d[1][j];
	g[2][j]+=d[1][j]*d[1][j];
      }
    }
    for(j=0;j<6;j++)
      for(k=0;k<3;k++)
	f[k][j]/=(double)group;
    
    for(j=0;j<6;j++)
      h[iter[0]][j]=(f[0][j]-(f[2][j]-f[1][j]*f[1][j])/T)/vol;
    
    mu[iter[0]]=2.0/9.0*(h[iter[0]][0]+h[iter[0]][1]+h[iter[0]][2])+h[iter[0]][3]+h[iter[0]][4]+h[iter[0]][5];
    mu[iter[0]]/=5.0;
    
    printf("\nf1, f2, f3 = %e %e %e\n",h[iter[0]][0],h[iter[0]][1],h[iter[0]][2]);
    printf("f4, f5, f6 = %e %e %e\n",h[iter[0]][3],h[iter[0]][4],h[iter[0]][5]);
    printf("group, mu = %d %e\n",iter[0],mu[iter[0]]);
    
  }
  
  sprintf(filename,"%s",argv[1]);
  kk=strlen(filename);
  close_xv4a_(filename,&PART,&time4,&pos[0],&vel[0],kk);

  for(j=0;j<6;j++)
    for(k=0;k<3;k++)
      g[k][j]/=(4.0*(double)group);
  
  for(j=0;j<6;j++)
    h[4][j]=(g[0][j]-(g[2][j]-g[1][j]*g[1][j])/T)/vol;

  mu[4]=2.0/9.0*(h[4][0]+h[4][1]+h[4][2])+h[4][3]+h[4][4]+h[4][5];
  mu[4]/=5.0;
  
  printf("\nf1, f2, f3 = %e %e %e\n",h[4][0],h[4][1],h[4][2]);
  printf("f4, f5, f6 = %e %e %e\n",h[4][3],h[4][4],h[4][5]);
  printf("final mu = %e\n",mu[4]);
  
}
 

double comppot(int class, double dx, int axis, float* pos, double* charge) {
  int j, k, n, x, y, z;
  double stretch[2], s[3], r, xx, yy, zz, potential;
  float tpos[PART][3];

  potential=0.0;

  if(class==0) {
    stretch[0]=1.0+dx+3.0/4.0*dx*dx;
    stretch[1]=1.0-0.5*dx;
    for(j=0;j<3;j++)
      s[j]=stretch[1];
    s[axis]=stretch[0];
    for(j=0;j<PART;j++) {
      for(n=0;n<3;n++)
	tpos[j][n]=pos[3*j+n]*s[n];
    }
    for(j=0;j<PART-1;j++){
      for(k=j+1;k<PART;k++) {
	for(x=0;x<3;x++) {
	  for(y=0;y<3;y++) {
	    for(z=0;z<3;z++) {
	      xx=tpos[j][0]-tpos[k][0]+(x-1)*BOXSIZE*s[0];
	      yy=tpos[j][1]-tpos[k][2]+(y-1)*BOXSIZE*s[1];
	      zz=tpos[j][2]-tpos[k][1]+(z-1)*BOXSIZE*s[2];
	      r=sqrt(xx*xx+yy*yy+zz*zz);
	      potential+=AHC*charge[j]*charge[k]*exp(-r/lambda)/r;
	    }
	  }
	}
      }
    }
  }
  else {

  }
  
}
