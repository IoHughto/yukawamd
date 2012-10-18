#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
  int j, k, l, num;
  double stretch, pos[3];
  FILE *out;

  out=fopen("bcc.config","w");

  num=atoi(argv[1]);

  stretch=38.1925;

  for(j=0;j<num;j++) {
    for(k=0;k<num;k++) {
      for(l=0;l<num;l++) {
	pos[0]=stretch*(double)j;
	pos[1]=stretch*(double)k;
	pos[2]=stretch*(double)l;
	fwrite(pos,sizeof(double),sizeof(pos)/sizeof(double),out);
      }
    }
  }  
  for(j=0;j<num;j++) {
    for(k=0;k<num;k++) {
      for(l=0;l<num;l++) {
	pos[0]=stretch*(double)j+stretch/2.0;
	pos[1]=stretch*(double)k+stretch/2.0;
	pos[2]=stretch*(double)l;
	fwrite(pos,sizeof(double),sizeof(pos)/sizeof(double),out);
      }
    }
  }  
  for(j=0;j<num;j++) {
    for(k=0;k<num;k++) {
      for(l=0;l<num;l++) {
	pos[0]=stretch*(double)j+stretch/2.0;
	pos[1]=stretch*(double)k;
	pos[2]=stretch*(double)l+stretch/2.0;
	fwrite(pos,sizeof(double),sizeof(pos)/sizeof(double),out);
      }
    }
  }  
  for(j=0;j<num;j++) {
    for(k=0;k<num;k++) {
      for(l=0;l<num;l++) {
	pos[0]=stretch*(double)j;
	pos[1]=stretch*(double)k+stretch/2.0;
	pos[2]=stretch*(double)l+stretch/2.0;
	fwrite(pos,sizeof(double),sizeof(pos)/sizeof(double),out);
      }
    }
  }
  fclose(out);
}
