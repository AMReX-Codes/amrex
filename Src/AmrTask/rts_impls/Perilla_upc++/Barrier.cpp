#include <Barrier.H>
#include <stdio.h>
#include <limits.h>
#include<assert.h>

Barrier::Barrier()
{
  //With this intializer, numthreads has to be specified when syncing, i.e. sync(numthreads)
  counter = INT_MAX;
  maxThreads=INT_MAX;
  condition= PTHREAD_COND_INITIALIZER;
  condition_mutex= PTHREAD_MUTEX_INITIALIZER;
}

Barrier::Barrier(int numthreads)
{
//With this initializer, both sync() and sync(numthreads) can be used
  counter = numthreads;
  maxThreads= numthreads;
  condition= PTHREAD_COND_INITIALIZER;
  condition_mutex= PTHREAD_MUTEX_INITIALIZER;
}

void Barrier::init(int numthreads)
{
//Similar to Barrier(int numthreads)
  counter = numthreads;
  maxThreads= numthreads;
  condition= PTHREAD_COND_INITIALIZER;
  condition_mutex= PTHREAD_MUTEX_INITIALIZER;
}

void Barrier::sync() //sync all threads associated with this barrier
{
  assert(maxThreads<INT_MAX);

  pthread_mutex_lock(&condition_mutex);
  counter--;

  if(counter == 0)
  {
      counter = maxThreads;
      pthread_cond_broadcast(&condition);
  } else {
      pthread_cond_wait(&condition, &condition_mutex);
  }
  pthread_mutex_unlock(&condition_mutex);
}

void Barrier::sync(int numthreads) //sync a subset of threads
{
  assert(numthreads<=maxThreads);
  pthread_mutex_lock(&condition_mutex);
  counter--;

  if(counter == (maxThreads-numthreads))
  {
      counter = maxThreads;
      pthread_cond_broadcast(&condition);
  } else {
      pthread_cond_wait(&condition, &condition_mutex);
  }
  pthread_mutex_unlock(&condition_mutex);
}
