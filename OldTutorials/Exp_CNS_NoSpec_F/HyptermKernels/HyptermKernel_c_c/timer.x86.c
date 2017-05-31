#include <stdint.h>

//inline uint64_t CycleTime(){
uint64_t CycleTime(){
  uint64_t lo, hi;
  __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
  return( (((uint64_t)hi) << 32) | ((uint64_t)lo) );
}
