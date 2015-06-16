#ifndef _TIMER_H_
#define _TIMER_H_
#include <stdio.h>
#include <time.h>
#include <sys/time.h>

inline
double dgettime(void) { // return current time in second
    struct timeval tv;
    struct timezone tz;
    gettimeofday(&tv, &tz);
    return (double)tv.tv_sec + (double)tv.tv_usec/ 1000000.0;
}

#endif
