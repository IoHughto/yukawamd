#include <stdio.h>

int main(int argc, char *argv[])
{
  int PART, n, skip, conf, j;
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
  
  double pos[PART][3], vel[PART][3];
  
  fprintf(out,"%d\n[name]\n",PART);
  conf=0;
  while(conf<skip)
    {
      for(j=0;j<PART;j++) {
	fread(pos[j],sizeof(double),3,in);
	fread(vel[j],sizeof(double),3,in);
      }
      conf++;
    }
  for(j=0;j<PART;j++) {
    fread(pos[j],sizeof(double),3,in);
    fread(vel[j],sizeof(double),3,in);
    fprintf(out,"C %lf %lf %lf\n",pos[j][0],pos[j][1],pos[j][2]);
  }
}
