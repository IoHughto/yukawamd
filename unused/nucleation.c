#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

void read_xv8_(char*,int*,double*,double*,double*,int);
void write_xv8_(char*,int*,double*,double*,double*,int);

int main(int argc, char *argv[]) {
  int j, k, l, n, PART, qflag, flag, difference;
  double rcut, rho, BOXSIZE, time8, diff[3], r, delta;
  FILE *inion, *oution;
  char filename[100], line[300];

  if(argc != 8) {
    fprintf(stderr,"Syntax Error: %s <ckpt in> <ion.dat in> <density> <r> <ckpt out> <ion.dat out> <quit if no match?>\n",argv[0]);
    return 1;
  }

  inion=fopen(argv[2],"r");
  sscanf(argv[3],"%lf",&rho);
  sscanf(argv[4],"%lf",&rcut);
  oution=fopen(argv[6],"w");
  qflag=atoi(argv[7]);

  fgets(line,300,inion);
  sscanf(line,"%d",&PART);
  
  int cutparticles[PART], newparticles[PART], ncut, nnew;
  double pos[PART][3], vel[PART][3], num;
  double crystal[PART][3];
  double charge[PART], mass[PART];

  BOXSIZE=pow(PART/rho,1.0/3.0);
  delta=pow(2.0/rho,1.0/3.0);
  num=pow(PART/2.0,1.0/3.0);

  j=0;

  while(fgets(line,300,inion)!=NULL) {
    sscanf(line,"%d %lf %lf\n",&k,&charge[j],&mass[j]);
    j++;
  }

  sprintf(filename,"%s",argv[1]);
  k=strlen(filename);
  read_xv8_(filename,&PART,&time8,&pos[0][0],&vel[0][0],k);
  
  ncut=0;
  for(j=0;j<PART;j++) {
    for(n=0;n<3;n++) {
      diff[n]=BOXSIZE/2.0-pos[j][n];
    }
    r=sqrt(diff[0]*diff[0]+diff[1]*diff[1]+diff[2]*diff[2]);
    if(r<rcut) {
      cutparticles[ncut]=j;
      ncut++;
    }
  }

  n=0;
  for(j=0;j<num;j++) {
    for(k=0;k<num;k++) {
      for(l=0;l<num;l++) {
	crystal[n][0]=delta*(double)j;
	crystal[n][1]=delta*(double)k;
	crystal[n][2]=delta*(double)l;
	n++;
      }
    }
  }  
  for(j=0;j<num;j++) {
    for(k=0;k<num;k++) {
      for(l=0;l<num;l++) {
	crystal[n][0]=delta*(double)j+delta/2.0;
	crystal[n][1]=delta*(double)k+delta/2.0;
	crystal[n][2]=delta*(double)l+delta/2.0;
	n++;
      }
    }
  }

  nnew=0;
  for(j=0;j<PART;j++) {
    for(n=0;n<3;n++) {
      diff[n]=BOXSIZE/2.0-crystal[j][n]+0.01;
    }
    r=sqrt(diff[0]*diff[0]+diff[1]*diff[1]+diff[2]*diff[2]);
    if(r<rcut) {
      newparticles[nnew]=j;
      nnew++;
    }
  }

  difference=abs(nnew-ncut);
  if(difference==0)
    difference++;

  int dparticles[difference];
  double rmax;

  if(nnew==ncut) {
    for(j=0;j<ncut;j++) {
      for(n=0;n<3;n++) {
	pos[cutparticles[j]][n]=crystal[newparticles[j]][n];
	vel[cutparticles[j]][n]*=sqrt(16.0/1.0e+24);
	mass[cutparticles[j]]=1.0e+24;
      }
    }
  }
  else if(qflag==1) {
    printf("Size mismatch!\n\nCut particles: %d\nCrystal particles: %d\n",ncut,nnew);
    return 1;
  }
  else {
    if(ncut>nnew) {
      for(l=0;l<difference;l++) {
	rmax=0.0;
	for(j=0;j<ncut;j++) {
	  for(n=0;n<3;n++) {
	    diff[n]=BOXSIZE/2.0-pos[cutparticles[j]][n]+0.01;
	  }
	  r=sqrt(diff[0]*diff[0]+diff[1]*diff[1]+diff[2]*diff[2]);
	  flag=0;
	  for(k=l-1;k>=0;k--) {
	    if(cutparticles[j]==dparticles[k])
	      flag=1;
	  }
	  if(flag==0) {
	    if(r>rmax) {
	      rmax=r;
	      dparticles[l]=cutparticles[j];
	    }
	  }
	}
      }
      for(j=0;j<nnew;j++) {
	flag=0;
	for(l=0;l<difference;l++) {
	  if(cutparticles[j]==dparticles[l])
	    flag==1;
	}
	if(flag==0) {
	  for(n=0;n<3;n++) {
	    pos[cutparticles[j]][n]=crystal[newparticles[j]][n];
	    vel[cutparticles[j]][n]*=sqrt(16.0/1.0e+24);
	    mass[cutparticles[j]]=1.0e+24;
	  }
	}
      }
    }
    else {
      for(l=0;l<difference;l++) {
	rmax=10000.0;
	for(j=0;j<PART;j++) {
	  flag=0;
	  for(k=0;k<ncut;k++) {
	    if(j==cutparticles[k])
	      flag=1;
	  }
	  for(k=l-1;k>=0;k--) {
	    if(j==dparticles[k])
	      flag=1;
	  }
	  if(flag==0) {
	    for(n=0;n<3;n++) {
	      diff[n]=BOXSIZE/2.0-pos[j][n]+0.01;
	    }
	    r=sqrt(diff[0]*diff[0]+diff[1]*diff[1]+diff[2]*diff[2]);
	    if(flag==0) {
	      if(r<rmax) {
		rmax=r;
		dparticles[l]=j;
	      }
	    }
	  }
	}
      }
      for(j=0;j<ncut;j++) {
	for(n=0;n<3;n++) {
	  pos[cutparticles[j]][n]=crystal[newparticles[j]][n];
	  vel[cutparticles[j]][n]*=sqrt(16.0/1.0e+24);
	  mass[cutparticles[j]]=1.0e+24;
	}
      }
      for(j=0;j<difference;j++){
	for(n=0;n<3;n++) {
	  pos[dparticles[j]][n]=crystal[newparticles[j+ncut]][n];
	  vel[dparticles[j]][n]*=sqrt(16.0/1.0e+24);
	  mass[dparticles[j]]=1.0e+24;
	}
      }
    }
  }

  sprintf(filename,"%s",argv[5]);
  k=strlen(filename);
  write_xv8_(filename,&PART,&time8,&pos[0][0],&vel[0][0],k);
  
  fprintf(oution,"%d\n",PART);
  for(j=0;j<PART;j++) {
    fprintf(oution,"%d %lf %lf\n",j,charge[j],mass[j]);
  }
}
