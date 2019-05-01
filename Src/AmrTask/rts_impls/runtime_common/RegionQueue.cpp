#include <RegionQueue.H>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

//////////////////////// class RegionQueue Definition Start /////////////////////////////////////  
RegionQueue::RegionQueue(void)
{
    max_size= perilla::TASK_QUEUE_DEFAULT_SIZE;
    buffer = new int[max_size];
    n = 0;
    front = 0;
    rear = 0;
    queueLock=PTHREAD_MUTEX_INITIALIZER;
}

RegionQueue::RegionQueue(int numTasks)
{
    buffer = new int[numTasks];
    n = 0;
    max_size = numTasks;
    front = 0;
    rear = 0;
    queueLock=PTHREAD_MUTEX_INITIALIZER;
}

RegionQueue::~RegionQueue()
{
    delete[] buffer;
}


void RegionQueue::addRegion(int r)
{
    pthread_mutex_lock(&queueLock);
    buffer[rear] = r;
    rear = (rear+1)%max_size;
    n++;
    pthread_mutex_unlock(&queueLock);
}

void RegionQueue::addRegion(int r, bool lockIgnored)
{
    if(!lockIgnored)pthread_mutex_lock(&queueLock);
    buffer[rear] = r;
    rear = (rear+1)%max_size;
    n++;
    if(!lockIgnored)pthread_mutex_unlock(&queueLock);
}

int RegionQueue::removeRegion()
{
    int r;
    pthread_mutex_lock(&queueLock);
    r = buffer[front];
    front = (front+1)%max_size;
    n--;
    pthread_mutex_unlock(&queueLock);
    return r;
}

int RegionQueue::removeRegion(bool lockIgnored)
{
    int r;
    if(!lockIgnored)pthread_mutex_lock(&queueLock);
    r = buffer[front];
    front = (front+1)%max_size;
    n--;
    if(!lockIgnored)pthread_mutex_unlock(&queueLock);
    return r;
}

int RegionQueue::getFrontRegion()
{
    return buffer[front];
}

int RegionQueue::getFrontRegion(bool lockIgnored)
{
    if(!lockIgnored)pthread_mutex_lock(&queueLock);
    return buffer[front];
    if(!lockIgnored)pthread_mutex_unlock(&queueLock);
}

int RegionQueue::queueSize()
{
    int size;
    pthread_mutex_lock(&queueLock);
    size = n;
    pthread_mutex_unlock(&queueLock);
    return size;
}

int RegionQueue::queueSize(bool lockIgnored)
{
    int size;
    if(!lockIgnored)pthread_mutex_lock(&queueLock);
    size = n;
    if(!lockIgnored)pthread_mutex_unlock(&queueLock);
    return size;
}
//////////////////////// class RegionQueue Definition End /////////////////////////////////////  
