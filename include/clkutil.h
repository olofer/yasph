/*
EXAMPLE:

  struct timespec tic, toc;

  clock_gettime(CLOCK_MONOTONIC, &tic);
  // do stuff
  clock_gettime(CLOCK_MONOTONIC, &toc);
  double elapsed_nanosecs = (double) timespecDiff(&tic, &toc);

*/

#ifndef __CLKUTIL_H__
#define __CLKUTIL_H__

int64_t timespecDiff(const struct timespec *timeA_p,
                     const struct timespec *timeB_p)
{
  return ((timeB_p->tv_sec * 1000000000) + timeB_p->tv_nsec) -
           ((timeA_p->tv_sec * 1000000000) + timeA_p->tv_nsec);
}

#define clkutil_stamp(C) clock_gettime(CLOCK_MONOTONIC, &(C))
#define clkutil_elapsed(A,B) timespecDiff(&(A), &(B))

#endif
