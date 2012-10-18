#include <stdio.h>

int main(int argc, char *argv[])
{
  int numskip, partnum, n, j;
  FILE *in, *out;

  if (argc != 5)
    {
      fprintf(stderr,"Syntax Error: %s <input> <output> <# of particles> <# of config to skip>\n",argv[0]);
      return 0;
    }
  
  numskip=atoi(argv[4]);
  partnum=atoi(argv[3]);
  
  double pos[3*partnum],vel[3*partnum];

  in=fopen(argv[1],"r");
  out=fopen(argv[2],"w");

  for(n=0;n<numskip;n++)
    {
      fread(pos,sizeof(double),3*partnum,in);
      fread(vel,sizeof(double),3*partnum,in);
    }
  fread(pos,sizeof(double),3*partnum,in);
  fread(vel,sizeof(double),3*partnum,in);
  fwrite(pos,sizeof(double),3*partnum,out);
  fwrite(vel,sizeof(double),3*partnum,out);
}
