#include  <stdlib.h>
#include  <unistd.h>
#include  <stdio.h>
#include  <string.h>
#include  <time.h>
#include  <sys/time.h>
#include  <sys/times.h>
#include  <sys/types.h>

/* Wallclock timer functions
 *   0 = reset
 *   1 = return seconds since last reset
 *   2 = return seconds since last call
 */

double          sec;    /* return value */
double          x;      /* current time */
double          x0;     /* time of previous call */
struct timeval  tv;
struct tms      cpu;
long            t0 = 0; /* helps retain precision of return value */

#pragma omp threadprivate(sec,x,x0,tv,cpu,t0)

double dwalltime00_(int *function)
{
  gettimeofday(&tv,0);
  if(*function==2) {
    x   = (tv.tv_sec - t0)  + tv.tv_usec/1000000.0;
    sec = x - x0;
  }
  else if(*function==1) {
    x   = (tv.tv_sec - t0)  + tv.tv_usec/1000000.0;
    sec = x;
  }
  else {
    t0 = tv.tv_sec;
    sec = x = 0.0;
  }
  x0 = x;
  return sec;
}
