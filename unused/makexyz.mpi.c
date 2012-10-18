#include <stdio.h>

int main(int argc, char *argv[])
{
  int PART, n, skip, conf;
  char line[300], inname[100], outname[100];
  FILE *in, *out;
  
  if(argc != 5)
    {
      fprintf(stderr,"Syntax Error: %s <input> <output> <# of particles> <# of configs to skip>\n",argv[0]);
      return 0;
    }
  
  in=fopen(argv[1],"r");
  out=fopen(argv[2],"w");
  
  skip=atoi(argv[4]);
  PART=atoi(argv[3]);
  
  double pos[3*PART], vel[3*PART];
  
  fprintf(out,"%d\n[name]\n",PART);
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
      fprintf(out,"C %lf %lf %lf\n",pos[3*n],pos[3*n+1],pos[3*n+2]);
    }
}
