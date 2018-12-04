
#include <RegionQueue.H>

//////////////////////// class RegionQueue Definition Start /////////////////////////////////////  
  RegionQueue::RegionQueue(void)
  {
    buffer = new int[perilla::TASK_QUEUE_DEFAULT_MAXSIZE];
    n = 0;
    bufSize = perilla::TASK_QUEUE_DEFAULT_MAXSIZE;
    front = 0;
    rear = 0;
    omp_init_lock(&queueLock);
  }

  RegionQueue::RegionQueue(int numTasks)
  {
    buffer = new int[numTasks];
    n = 0;
    bufSize = numTasks;
    front = 0;
    rear = 0;
    omp_init_lock(&queueLock);
  }

  void RegionQueue::addRegion(int r)
  {
    omp_set_lock(&queueLock);
    buffer[rear] = r;
    rear = (rear+1)%bufSize;
    n++;
    omp_unset_lock(&queueLock);
  }

  void RegionQueue::addRegion(int r, bool lockIgnore)
  {
    if(!lockIgnore)omp_set_lock(&queueLock);
    buffer[rear] = r;
    rear = (rear+1)%bufSize;
    n++;
    if(!lockIgnore)omp_unset_lock(&queueLock);
  }

  int RegionQueue::removeRegion()
  {
    int r;
    omp_set_lock(&queueLock);
    r = buffer[front];
    front = (front+1)%bufSize;
    n--;
    omp_unset_lock(&queueLock);
    return r;
  }

  int RegionQueue::removeRegion(bool lockIgnore)
  {
    int r;
    if(!lockIgnore)omp_set_lock(&queueLock);
    r = buffer[front];
    front = (front+1)%bufSize;
    n--;
    if(!lockIgnore)omp_unset_lock(&queueLock);
    return r;
  }

  int RegionQueue::getFrontRegion()
  {
    return buffer[front];
  }

  int RegionQueue::getFrontRegion(bool lockIgnore)
  {
    if(!lockIgnore)omp_set_lock(&queueLock);
    return buffer[front];
    if(!lockIgnore)omp_unset_lock(&queueLock);
  }

  int RegionQueue::queueSize()
  {
    int size;
    omp_set_lock(&queueLock);
    size = n;
    omp_unset_lock(&queueLock);
    return size;
  }

  int RegionQueue::queueSize(bool lockIgnore)
  {
    int size;
    if(!lockIgnore)omp_set_lock(&queueLock);
    size = n;
    if(!lockIgnore)omp_unset_lock(&queueLock);
    return size;
  }
//////////////////////// class RegionQueue Definition End /////////////////////////////////////  
