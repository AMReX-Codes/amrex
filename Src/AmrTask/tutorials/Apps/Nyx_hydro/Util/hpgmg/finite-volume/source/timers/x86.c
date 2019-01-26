//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
#include <stdint.h>
#define CALIBRATE_TIMER // mg.c will calibrate the timer to determine seconds per cycle
double getTime(){
  uint64_t lo, hi;
  __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
  return( 1e-9*((double)( (((uint64_t)hi) << 32) | ((uint64_t)lo) )) ); // timers are in units of seconds;  assume 1GHz cycle counter and convert later
}
