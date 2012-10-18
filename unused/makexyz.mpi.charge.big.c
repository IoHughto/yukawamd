#include <stdio.h>

int main(int argc, char *argv[])
{
  int PART, n, skip, conf;
  char line[300], inname[100], outname[100];
  FILE *in, *out, *chargefile;
  
  if(argc != 6)
    {
      fprintf(stderr,"Syntax Error: %s <positions> <charges> <output> <# of particles> <# of configs to skip>\n",argv[0]);
      return 0;
    }
  
  in=fopen(argv[1],"r");
  out=fopen(argv[3],"w");
  chargefile=fopen(argv[2],"r");
  
  skip=atoi(argv[5]);
  PART=atoi(argv[4]);
  
  double pos[3*PART], vel[3*PART], charge[PART], BOXSIZE, rho;

  BOXSIZE=363.761;

  for(n=0;n<PART;n++) {
    fgets(line,300,chargefile);
    sscanf(line,"%lf",&charge[n]);
  }

  fprintf(out,"%d\n[name]\n",PART*8);
  conf=0;
  while(conf<skip)
    {
      fread(pos,sizeof(double),3*PART,in);
      fread(vel,sizeof(double),3*PART,in);
      conf++;
    }
  fread(pos,sizeof(double),3*PART,in);
  fread(vel,sizeof(double),3*PART,in);
  for(n=0;n<PART;n++)
    {
      if(charge[n]<7.0) {
	fprintf(out,"C %lf %lf %lf\n",pos[3*n],pos[3*n+1],pos[3*n+2]);
	fprintf(out,"C %lf %lf %lf\n",pos[3*n]+BOXSIZE,pos[3*n+1],pos[3*n+2]);
	fprintf(out,"C %lf %lf %lf\n",pos[3*n],pos[3*n+1]+BOXSIZE,pos[3*n+2]);
	fprintf(out,"C %lf %lf %lf\n",pos[3*n],pos[3*n+1],pos[3*n+2]+BOXSIZE);
	fprintf(out,"C %lf %lf %lf\n",pos[3*n]+BOXSIZE,pos[3*n+1]+BOXSIZE,pos[3*n+2]);
	fprintf(out,"C %lf %lf %lf\n",pos[3*n]+BOXSIZE,pos[3*n+1],pos[3*n+2]+BOXSIZE);
	fprintf(out,"C %lf %lf %lf\n",pos[3*n],pos[3*n+1]+BOXSIZE,pos[3*n+2]+BOXSIZE);
	fprintf(out,"C %lf %lf %lf\n",pos[3*n]+BOXSIZE,pos[3*n+1]+BOXSIZE,pos[3*n+2]+BOXSIZE);
      }
      else {
	fprintf(out,"O %lf %lf %lf\n",pos[3*n],pos[3*n+1],pos[3*n+2]);
	fprintf(out,"O %lf %lf %lf\n",pos[3*n]+BOXSIZE,pos[3*n+1],pos[3*n+2]);
	fprintf(out,"O %lf %lf %lf\n",pos[3*n],pos[3*n+1]+BOXSIZE,pos[3*n+2]);
	fprintf(out,"O %lf %lf %lf\n",pos[3*n],pos[3*n+1],pos[3*n+2]+BOXSIZE);
	fprintf(out,"O %lf %lf %lf\n",pos[3*n]+BOXSIZE,pos[3*n+1]+BOXSIZE,pos[3*n+2]);
	fprintf(out,"O %lf %lf %lf\n",pos[3*n]+BOXSIZE,pos[3*n+1],pos[3*n+2]+BOXSIZE);
	fprintf(out,"O %lf %lf %lf\n",pos[3*n],pos[3*n+1]+BOXSIZE,pos[3*n+2]+BOXSIZE);
	fprintf(out,"O %lf %lf %lf\n",pos[3*n]+BOXSIZE,pos[3*n+1]+BOXSIZE,pos[3*n+2]+BOXSIZE);
      }
    }
}
