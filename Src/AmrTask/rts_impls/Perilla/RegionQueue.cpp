#include <RegionQueue.H>

//////////////////////// class RegionQueue Definition Start /////////////////////////////////////  
  RegionQueue::RegionQueue(void)
  {
    buffer = new int[perilla::TASK_QUEUE_DEFAULT_MAXSIZE];
    n = 0;
    bufSize = perilla::TASK_QUEUE_DEFAULT_MAXSIZE;
    front = 0;
    rear = 0;
    queueLock=PTHREAD_MUTEX_INITIALIZER;
  }

  RegionQueue::RegionQueue(int numTasks)
  {
    buffer = new int[numTasks];
    n = 0;
    bufSize = numTasks;
    front = 0;
    rear = 0;
    queueLock=PTHREAD_MUTEX_INITIALIZER;
  }

  void RegionQueue::addRegion(int r)
  {
    pthread_mutex_lock(&queueLock);
    buffer[rear] = r;
    rear = (rear+1)%bufSize;
    n++;
    pthread_mutex_unlock(&queueLock);
  }

  void RegionQueue::addRegion(int r, bool lockIgnore)
  {
    if(!lockIgnore)pthread_mutex_lock(&queueLock);
    buffer[rear] = r;
    rear = (rear+1)%bufSize;
    n++;
    if(!lockIgnore)pthread_mutex_unlock(&queueLock);
  }

  int RegionQueue::removeRegion()
  {
    int r;
    pthread_mutex_lock(&queueLock);
    r = buffer[front];
    front = (front+1)%bufSize;
    n--;
    pthread_mutex_unlock(&queueLock);
    return r;
  }

  int RegionQueue::removeRegion(bool lockIgnore)
  {
    int r;
    if(!lockIgnore)pthread_mutex_lock(&queueLock);
    r = buffer[front];
    front = (front+1)%bufSize;
    n--;
    if(!lockIgnore)pthread_mutex_unlock(&queueLock);
    return r;
  }

  int RegionQueue::getFrontRegion()
  {
    return buffer[front];
  }

  int RegionQueue::getFrontRegion(bool lockIgnore)
  {
    if(!lockIgnore)pthread_mutex_lock(&queueLock);
    return buffer[front];
    if(!lockIgnore)pthread_mutex_unlock(&queueLock);
  }

  int RegionQueue::queueSize()
  {
    int size;
    pthread_mutex_lock(&queueLock);
    size = n;
    pthread_mutex_unlock(&queueLock);
    return size;
  }

  int RegionQueue::queueSize(bool lockIgnore)
  {
    int size;
    if(!lockIgnore)pthread_mutex_lock(&queueLock);
    size = n;
    if(!lockIgnore)pthread_mutex_unlock(&queueLock);
    return size;
  }
//////////////////////// class RegionQueue Definition End /////////////////////////////////////  
