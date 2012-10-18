#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define min(x1,x2) ((x1) > (x2) ? (x2):(x1))

void read_xv8_(char*,int*,double*,double*,double*,int);
double findshort(double a, double b, double *lr, double BOXSIZE);

#define PI 3.141592653589
#define ALPHA 1.0/137.035999
#define HBARC 197.32693

int main(int argc, char *argv[]) {

  int PART, binnum, n, j, k, l, m, iter, flag, ionnum[5], dt, dumpi, ierr, rank, count, compnum, conf, slices, sflag;
  long startt, stopt, confnum;
  double rho, BOXSIZE, dump, relpos[3], r, bind, time8, lr[3], rp[3], sliced;
  char line[300], filename[100];
  FILE *out, *ionfile;

  if(argc != 10) {
      fprintf(stderr,"Syntax Error: %s <ion.dat file> <output base> <# of particles> <density> <# of slices> <# of bins> <start time> <end time> <data cadence>\n",argv[0]);
      return 1;
    }

  ionfile=fopen(argv[1],"r");
  PART=atoi(argv[3]);
  sscanf(argv[4],"%lf",&rho);
  slices=atoi(argv[5]);
  binnum=atoi(argv[6]);
  startt=atol(argv[7]);
  stopt=atol(argv[8]);
  dt=atoi(argv[9]);

  FILE *xyz[slices];
  for(j=0;j<slices;j++) {
    sprintf(filename,"nn.%d.xyz",j+1);
    xyz[j]=fopen(filename,"w");
  }

  confnum=(stopt-startt)/dt+1;

  BOXSIZE=pow(PART/rho,1.0/3.0);
  bind=PI/binnum;
  sliced=BOXSIZE/slices;

  fgets(line,300,ionfile);

  double charge[PART], mass[PART];

  for(j=0;j<5;j++)
    ionnum[j]=0;

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

  double vel[3*PART], pos[3*PART];
  double rmin[14], x[14], y[14], z[14], theta;
  int bins[slices][binnum], snum[slices];

  for(j=0;j<slices;j++)
    snum[j]=0;

  for(l=0;l<binnum;l++) 
    for(j=0;j<slices;j++)
      bins[j][l]=0;
  
  for(conf=0;conf<confnum;conf++) {
    sprintf(filename,"md.ckpt.%011d",startt+conf*dt);
    k=strlen(filename);
    printf("%s\n",filename);
    read_xv8_(filename,&PART,&time8,&pos[0],&vel[0],k);
    for(j=0;j<PART;j++) {
      for(n=0;n<14;n++)
	rmin[n]=100.0;
      for(k=0;k<PART;k++) {
	if(k!=j) {
	  for(n=0;n<3;n++) {
	    rp[n]=findshort(pos[3*j+n],pos[3*k+n],&lr[n],BOXSIZE);
	    rp[n]*=lr[n];
	  }
	  r=sqrt(rp[0]*rp[0]+rp[1]*rp[1]+rp[2]*rp[2]);
	  if(r<rmin[13]) {
	    if(r<rmin[12]) {
	      if(r<rmin[11]) {
		if(r<rmin[10]) {
		  if(r<rmin[9]) {
		    if(r<rmin[8]) {
		      if(r<rmin[7]) {
			if(r<rmin[6]) {
			  if(r<rmin[5]) {
			    if(r<rmin[4]) {
			      if(r<rmin[3]) {
				if(r<rmin[2]) {
				  if(r<rmin[1]) {
				    if(r<rmin[0]) {
				      for(n=13;n>0;n--) {
					rmin[n]=rmin[n-1];
					x[n]=x[n-1];
					y[n]=y[n-1];
					z[n]=z[n-1];
				      }
				      rmin[0]=r;
				      x[0]=rp[0];
				      y[0]=rp[1];
				      z[0]=rp[2];
				    }
				    else {
				      for(n=13;n>1;n--) {
					rmin[n]=rmin[n-1];
					x[n]=x[n-1];
					y[n]=y[n-1];
					z[n]=z[n-1];
				      }
				      rmin[1]=r;
				      x[1]=rp[0];
				      y[1]=rp[1];
				      z[1]=rp[2];
				    }
				  }
				  else {
				    for(n=13;n>2;n--) {
				      rmin[n]=rmin[n-1];
				      x[n]=x[n-1];
				      y[n]=y[n-1];
				      z[n]=z[n-1];
				    }
				    rmin[2]=r;
				    x[2]=rp[0];
				    y[2]=rp[1];
				    z[2]=rp[2];
				  }
				}
				else {
				  for(n=13;n>3;n--) {
				    rmin[n]=rmin[n-1];
				    x[n]=x[n-1];
				    y[n]=y[n-1];
				    z[n]=z[n-1];
				  }
				  rmin[3]=r;
				  x[3]=rp[0];
				  y[3]=rp[1];
				  z[3]=rp[2];
				}
			      }
			      else {
				for(n=13;n>4;n--) {
				  rmin[n]=rmin[n-1];
				  x[n]=x[n-1];
				  y[n]=y[n-1];
				  z[n]=z[n-1];
				}
				rmin[4]=r;
				x[4]=rp[0];
				y[4]=rp[1];
				z[4]=rp[2];
			      }			
			    }
			    else {
			      for(n=13;n>5;n--) {
				rmin[n]=rmin[n-1];
				x[n]=x[n-1];
				y[n]=y[n-1];
				z[n]=z[n-1];
			      }
			      rmin[5]=r;
			      x[5]=rp[0];
			      y[5]=rp[1];
			      z[5]=rp[2];
			    }
			  }
			  else {
			    for(n=13;n>6;n--) {
			      rmin[n]=rmin[n-1];
			      x[n]=x[n-1];
			      y[n]=y[n-1];
			      z[n]=z[n-1];
			    }
			    rmin[6]=r;
			    x[6]=rp[0];
			    y[6]=rp[1];
			    z[6]=rp[2];
			  }
			}
			else {
			  for(n=13;n>7;n--) {
			    rmin[n]=rmin[n-1];
			    x[n]=x[n-1];
			    y[n]=y[n-1];
			    z[n]=z[n-1];
			  }
			  rmin[7]=r;
			  x[7]=rp[0];
			  y[7]=rp[1];
			  z[7]=rp[2];
			}
		      }
		      else {
			for(n=13;n>8;n--) {
			  rmin[n]=rmin[n-1];
			  x[n]=x[n-1];
			  y[n]=y[n-1];
			  z[n]=z[n-1];
			}
			rmin[8]=r;
			x[8]=rp[0];
			y[8]=rp[1];
			z[8]=rp[2];
		      }
		    }
		    else {
		      for(n=13;n>9;n--) {
			rmin[n]=rmin[n-1];
			x[n]=x[n-1];
			y[n]=y[n-1];
			z[n]=z[n-1];
		      }
		      rmin[9]=r;
		      x[9]=rp[0];
		      y[9]=rp[1];
		      z[9]=rp[2];
		    }
		  }
		  else {
		    for(n=13;n>10;n--) {
		      rmin[n]=rmin[n-1];
		      x[n]=x[n-1];
		      y[n]=y[n-1];
		      z[n]=z[n-1];
		    }
		    rmin[10]=r;
		    x[10]=rp[0];
		    y[10]=rp[1];
		    z[10]=rp[2];
		  }
		}
		else {
		  for(n=13;n>11;n--) {
		    rmin[n]=rmin[n-1];
		    x[n]=x[n-1];
		    y[n]=y[n-1];
		    z[n]=z[n-1];
		  }
		  rmin[11]=r;
		  x[11]=rp[0];
		  y[11]=rp[1];
		  z[11]=rp[2];
		}
	      }
	      else {
		for(n=13;n>12;n--) {
		  rmin[n]=rmin[n-1];
		  x[n]=x[n-1];
		  y[n]=y[n-1];
		  z[n]=z[n-1];
		}
		rmin[12]=r;
		x[12]=rp[0];
		y[12]=rp[1];
		z[12]=rp[2];
	      }
	    }
	    else {
	      rmin[13]=r;
	      x[13]=rp[0];
	      y[13]=rp[1];
	      z[13]=rp[2];
	    }
	  }
	}
      }
      sflag=0;
      for(m=0;m<slices;m++) {
	if(sflag==0) {
	  if(pos[3*j+2]<sliced*(m+1)) {
	    snum[m]++;
	    sflag=1;
	    for(k=0;k<14;k++) {
	      fprintf(xyz[m],"C %e %e %e\n",x[k],y[k],z[k]);
	      for(l=0;l<14;l++) {
		if(k!=l) {
		  theta=x[k]*x[l]+y[k]*y[l]+z[k]*z[l];
		  theta/=rmin[k];
		  theta/=rmin[l];
		  theta=acos(theta);
		  flag=0;
		  for(n=0;n<binnum;n++) {
		    if(flag==0) {
		      if(theta<bind*(n+1)) {
			bins[m][n]++;
			flag=1;
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  
  sprintf(filename,"%s.nn.out",argv[2]);
  out=fopen(filename,"w");
  for(n=0;n<binnum;n++) {
    fprintf(out,"%e ",n*bind);
    for(m=0;m<slices;m++) {
      fprintf(out,"%e ",(double)bins[m][n]/(56.0*snum[m]*confnum));
    }
    fprintf(out,"\n");
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
