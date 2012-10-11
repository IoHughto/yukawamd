#include  <stdlib.h>
#include  <unistd.h>
#include  <stdio.h>
#include  <string.h>
#include  <time.h>
#include  <sys/time.h>
#include  <sys/times.h>
#include  <sys/types.h>

/*
 * Arguments:
 *   double *dt   (in/out) -- the total accumulated by the timer
 *   double *self (in/out) -- total minus time accumulated by enclosed timers
 *   int    *n    (in/out) -- counter
 *
 * Arguments dt, self and n must be maintained by the user. Results are accumulated
 * into these variables.
 */

typedef struct {
  double  time;
  double  others;
} xtimer;
xtimer  xclock[21];     /* timer stack */

xtimer          *p;

struct timeval  tv;
struct tms      cpu;
long            t0;     /* helps retain precision of return value */

#pragma omp threadprivate(p,xclock,tv,cpu,t0)


/***********************************************************************
 *  Call inittime once at the beginning of a run, to initialize a pointer to
 *  the top of the timer stack, and to record an initial time, t0. Time t0
 *  will be subtracted from all succeeding measurements.
 */
void inittime_(void)
{
  p = xclock;
  gettimeofday(&tv,0);
  t0 = tv.tv_sec;
}


/***********************************************************************
 *  Start a timer, and push it onto the timer stack.
 */
void starttimer_(void)
{
  p++;
  gettimeofday(&tv,0);
  p->time = (tv.tv_sec-t0) + tv.tv_usec/1000000.0;
  p->others = 0.0;
}


/***********************************************************************
 *  Stop the current timer on the top of the stack, and pop it.
 */
void stoptimer_(int *accumulate, double *deltat, double *self, int *n)
{
  double   x;
  gettimeofday(&tv,0);
  p->time = ( (tv.tv_sec-t0) + tv.tv_usec/1000000.0 ) - p->time;
  x = p->time;
  if(accumulate) {
    *deltat  += x;
    *self += x - p->others;
    *n += 1;
  }
  else {
    *deltat  = x;
    *self = x - p->others;
    *n = 1;
  }
  p--;
  p->others += x;
}
