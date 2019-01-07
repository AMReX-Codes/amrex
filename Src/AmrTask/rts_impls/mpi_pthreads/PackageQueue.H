#ifndef P_PACKAGEQUEUE_H
#define P_PACKAGEQUEUE_H

#include <PerillaConfig.H>
#include <pthread.h>
#include <mpi.h>

class Package
{
private:
  int source, destination;
public:
  double *databuf;
  int bufSize;
  pthread_mutex_t packageLock;
  volatile bool completed; //message transfer is done
  volatile bool served; //message transfer request has been served but may have not completed
  volatile bool notified;
  MPI_Request request; //!for MPI
  Package();
  ~Package();
  Package(int size);
  Package(int src, int dest);
  Package(int src, int dest, int size);
  void setPackageSource(int src);
  void setPackageDestination(int dest);
  void completeRequest(void);
  void completeRequest(bool lockIgnore);
  bool checkRequest(void);
  void generatePackage(int size);
};

class PackageQueue
{
private:
  Package *buffer[perilla::MSG_QUEUE_DEFAULT_MAXSIZE];
  volatile int n;
  volatile int front;
  volatile int rear;
  volatile int prear;
  int max_size;
public:  
  pthread_mutex_t queueLock;
  PackageQueue();
  ~PackageQueue();
  int queueSize(void);  
  int queueSize(bool lockIgnore);
  void enqueue(Package* package);  
  void enqueue(Package* package, bool lockIgnore);
  Package* dequeue(void);  
  Package* dequeue(bool lockIgnore);
  Package* getRear(void);
  Package* getRear(bool lockIgnore);
  Package* getFront(void);
  Package* getFront(bool lockIgnore);
  void emptyQueue();
};
#endif

