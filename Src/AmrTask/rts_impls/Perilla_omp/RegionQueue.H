#ifndef P_REGIONQUEUE_H
#define P_REGIONQUEUE_H

#include <PerillaConfig.H>
#include <omp.h>

//////////////////////// class RegionQueue Declaration Start /////////////////////////////////////
class RegionQueue
{
private:
  int* buffer;
  int bufSize;
  int n;
  int front;
  int rear;
  omp_lock_t queueLock;
public:
  RegionQueue();
  RegionQueue(int numTasks);
  void addRegion(int r);
  void addRegion(int r, bool lockIgnore);
  int removeRegion();
  int removeRegion(bool lockIgnore);  
  int getFrontRegion();
  int getFrontRegion(bool lockIgnore);
  int queueSize(bool lockIgnore);
  int queueSize();
};
//////////////////////// class RegionQueue Declaration End /////////////////////////////////////


#endif
