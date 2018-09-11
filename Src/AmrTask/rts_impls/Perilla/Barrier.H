#ifndef P_BARRIER_H
#define P_BARRIER_H
#include <pthread.h>

class Barrier
{
private:
  volatile int counter;
  int maxThreads;
  pthread_cond_t condition;
  pthread_mutex_t condition_mutex;
  volatile bool globalSense;
public:
  Barrier();
  Barrier(int maxThreads);
  void init(int maxThreads);
  void sync();
  void sync(int numthreads);
};
#endif
