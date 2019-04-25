#include <Barrier.H>
#include <stdio.h>
#include <limits.h>
#include<assert.h>

Barrier::Barrier()
{
  //With this intializer, numthreads has to be specified when syncing, i.e. sync(numthreads)
  counter = INT_MAX;
  globalSense = false;
  maxThreads=INT_MAX;
}

Barrier::Barrier(int numthreads)
{
//With this initializer, both sync() and sync(numthreads) can be used
#pragma omp critical
{
  counter = numthreads;
  maxThreads= numthreads;
  globalSense = false;
}
}

void Barrier::init(int numthreads)
{
//Similar to Barrier(int numthreads)
  counter = numthreads;
  maxThreads= numthreads;
  globalSense = false;
}

void Barrier::sync() //sync all threads associated with this barrier
{
  assert(maxThreads<INT_MAX);
  bool localSense;
  localSense = globalSense;
  localSense =  !localSense;
#pragma omp critical
  {
    counter--;
    if(counter == 0)
      {
        counter = maxThreads;
        globalSense = localSense;
      }
  }

  while(globalSense != localSense)
    {
#pragma omp flush
    }
}

void Barrier::sync(int numthreads) //sync a subset of threads
{
  assert(numthreads<=maxThreads);
  bool localSense;
  
  localSense = globalSense;  
  localSense =  !localSense;
  
#pragma omp critical
  {
    counter--;
    if(counter == (maxThreads-numthreads))
      {
	counter = maxThreads;
	globalSense = localSense;
      }
  }

  while(globalSense != localSense)
    {
#pragma omp flush
    }
}
