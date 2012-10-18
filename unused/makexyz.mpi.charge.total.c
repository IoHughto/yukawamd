#include <stdio.h>

int main(int argc, char *argv[])
{
  int PART, n, skip, conf;
  char line[300], inname[100], outname[100];
  FILE *in, *out, *chargefile;
  
  if(argc != 5)
    {
      fprintf(stderr,"Syntax Error: %s <positions> <charges> <# of particles> <# of configs>\n",argv[0]);
      return 0;
    }
  
  in=fopen(argv[1],"r");
  chargefile=fopen(argv[2],"r");
  
  skip=atoi(argv[4]);
  PART=atoi(argv[3]);
  
  double pos[3*PART], vel[3*PART], charge[PART];
  
  for(n=0;n<PART;n++) {
    fgets(line,300,chargefile);
    sscanf(line,"%lf",&charge[n]);
  }

  fclose(chargefile);

  conf=0;
  while(conf<skip)
    {
      fread(pos,sizeof(double),3*PART,in);
      fread(vel,sizeof(double),3*PART,in);
      sprintf(outname,"yukawa.%d.xyz",conf);
      out=fopen(outname,"w");
      fprintf(out,"%d\n[name]\n",PART);
      for(n=0;n<PART;n++) {
	if(charge[n]<7.0)
	  fprintf(out,"C %lf %lf %lf\n",pos[3*n],pos[3*n+1],pos[3*n+2]);
	else
	  fprintf(out,"O %lf %lf %lf\n",pos[3*n],pos[3*n+1],pos[3*n+2]);
      }
      fclose(out);
      conf++;
    }
  fclose(in);
}
