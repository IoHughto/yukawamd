#include  <stdlib.h>
#include  <stdio.h>
#include  <string.h>
#include  <sys/stat.h>
#include  <unistd.h>


void file_stat_(char *file, int *size, int k)
{
  struct stat   s;
  char   xfile[1025];
  (void)strncpy(xfile,file,k);
  xfile[k]=0;
/*printf("file_stat:  ***%s***\n", xfile);*/
  (void)stat(xfile,&s);
  *size = s.st_size;
  return;
}
