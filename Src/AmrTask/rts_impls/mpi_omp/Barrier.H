#ifndef P_BARRIER_H
#define P_BARRIER_H
#include <omp.h>

class Barrier
{
private:
  volatile int counter;
  int maxThreads;
  volatile bool globalSense;
public:
  Barrier();
  Barrier(int maxThreads);
  void init(int maxThreads);
  void sync();
  void sync(int numthreads);
};
#endif
