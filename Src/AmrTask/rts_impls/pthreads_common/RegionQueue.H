#ifndef P_REGIONQUEUE_H
#define P_REGIONQUEUE_H

#include <PerillaConfig.H>
#include <pthread.h>

//////////////////////// class RegionQueue Declaration Start /////////////////////////////////////
class RegionQueue
{
private:
  int* buffer;
  int n;
  int front;
  int rear;
  int max_size;
  pthread_mutex_t queueLock;
public:
  RegionQueue();
  RegionQueue(int numTasks);
  ~RegionQueue();
  void addRegion(int r);
  void addRegion(int r, bool lockIgnore);
  int removeRegion();
  int removeRegion(bool lockIgnore);  
  int getFrontRegion();
  int getFrontRegion(bool lockIgnore);
  int queueSize(bool lockIgnore);
  int queueSize();
  void reallocate();
};
//////////////////////// class RegionQueue Declaration End /////////////////////////////////////

#endif
