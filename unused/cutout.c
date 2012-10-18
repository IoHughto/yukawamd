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
  
  double pos[partnum][3],vel[partnum][3];

  in=fopen(argv[1],"r");
  out=fopen(argv[2],"w");

  for(n=0;n<numskip;n++)
    {
      for(j=0;j<partnum;j++)
	{
	  fread(pos[j],sizeof(double),3,in);
	  fread(vel[j],sizeof(double),3,in);
	}
    }
  for(j=0;j<partnum;j++)
    {
      fread(pos[j],sizeof(double),3,in);
      fread(vel[j],sizeof(double),3,in);
      fwrite(pos[j],sizeof(double),3,out);
      fwrite(vel[j],sizeof(double),3,out);
    }
}
