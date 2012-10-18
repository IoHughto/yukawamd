#include <stdio.h>

int main(int argc, char *argv[])
{
  int first, last, partnum, n, j;
  FILE *in, *out;

  if (argc != 6)
    {
      fprintf(stderr,"Syntax Error: %s <input> <output> <# of particles> <first config> <last config>\n",argv[0]);
      return 0;
    }
  
  partnum=atoi(argv[3]);
  last=atoi(argv[5]);
  first=atoi(argv[4]);
  
  double pos[3*partnum],vel[3*partnum];

  in=fopen(argv[1],"r");
  out=fopen(argv[2],"w");

  for(n=1;n<first;n++)
    {
      fread(pos,sizeof(double),3*partnum,in);
      fread(vel,sizeof(double),3*partnum,in);
    }
  for(n=first;n<=last;n++)
    {
      fread(pos,sizeof(double),3*partnum,in);
      fread(vel,sizeof(double),3*partnum,in);
      fwrite(pos,sizeof(double),3*partnum,out);
      fwrite(vel,sizeof(double),3*partnum,out);
    }
}
