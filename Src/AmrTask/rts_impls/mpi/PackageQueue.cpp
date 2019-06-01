#include <PackageQueue.H>
#include <iostream>
using namespace perilla;
#ifdef PERILLA_DEBUG
#include <PerillaMemCheck.H>
extern PerillaMemCheck memcheck;
#endif

Package::Package()
{
    databuf = 0;
    bufSize = 0;
    source = 0;
    destination = 0;
    completed = false;
    served = false;
    request = MPI_REQUEST_NULL;
    packageLock= PTHREAD_MUTEX_INITIALIZER;
#ifdef PERILLA_DEBUG
    memcheck.add(memcheck.genKey(this), (void*)this, "Package");
#endif
}

Package::~Package()
{
    if(databuf) free(databuf);
#ifdef PERILLA_DEBUG
    memcheck.remove(memcheck.genKey(this));
#endif
}

Package::Package(int size)
{
    databuf = new double[size];
    bufSize = size;
    source = 0;
    destination = 0;
    completed = false;
    served = false;
    request = MPI_REQUEST_NULL;
    packageLock= PTHREAD_MUTEX_INITIALIZER;
#ifdef PERILLA_DEBUG
    memcheck.add(memcheck.genKey(this), (void*)this, "Package");
#endif
}

Package::Package(int src, int dest)
{
    bufSize = 0;
    source = src;
    destination = dest;
    completed = false;
    served = false;
    request = MPI_REQUEST_NULL;
    packageLock= PTHREAD_MUTEX_INITIALIZER;
#ifdef PERILLA_DEBUG
    memcheck.add(memcheck.genKey(this), (void*)this, "Package");
#endif
}

Package::Package(int src, int dest, int size)
{
    databuf = new double[size];
    bufSize = size;
    source = src;
    destination = dest;
    completed = false;
    served = false;
    request = MPI_REQUEST_NULL;
    packageLock= PTHREAD_MUTEX_INITIALIZER;
#ifdef PERILLA_DEBUG
    memcheck.add(memcheck.genKey(this), (void*)this, "Package");
#endif
}

void Package::setPackageSource(int src)
{
    source = src;
}

void Package::setPackageDestination(int dest)
{
    destination = dest;
}

void Package::completeRequest(void)
{
    pthread_mutex_lock(&packageLock);
    completed = true;
    pthread_mutex_unlock(&packageLock);
}

void Package::completeRequest(bool canAvoidLock)
{
    if(!canAvoidLock)pthread_mutex_lock(&packageLock);
    completed = true;
    if(!canAvoidLock)pthread_mutex_unlock(&packageLock);
}

bool Package::checkRequest(void)
{
    return completed;
}  

void Package::generatePackage(int size)
{
    databuf = new double[size];
    bufSize = size;
    source = 0;
    destination = 0;
    completed = false;
    served = false;
    request = MPI_REQUEST_NULL;
    packageLock= PTHREAD_MUTEX_INITIALIZER;
#ifdef PERILLA_DEBUG
    memcheck.add(memcheck.genKey(this), (void*)this, "Package");
#endif
}

PackageQueue::PackageQueue()
{
    n = 0;
    front = 0;
    rear = 0;
    prear = -1;
    queueLock= PTHREAD_MUTEX_INITIALIZER;;
}

int PackageQueue::queueSize(void)
{
    int size;
    pthread_mutex_lock(&queueLock);
    size = n;
    pthread_mutex_unlock(&queueLock);
    return size;
}

int PackageQueue::queueSize(bool canAvoidLock)
{
    int size;
    if(!canAvoidLock)pthread_mutex_lock(&queueLock);
    size = n;
    if(!canAvoidLock)pthread_mutex_unlock(&queueLock);
    return size;
}

void PackageQueue::enqueue(Package* package)
{
    pthread_mutex_lock(&queueLock);
#ifdef PERILLA_DEBUG
    if(n==perilla::MSG_QUEUE_DEFAULT_MAXSIZE){
        printf("Failed to Enqueue: Queue Overflow\n");
        exit(0);
    }
#endif
    buffer[rear] = package;
    prear = rear;    
    rear = (rear+1)%perilla::MSG_QUEUE_DEFAULT_MAXSIZE;
    n++;
    pthread_mutex_unlock(&queueLock);
}

void PackageQueue::enqueue(Package* package, bool canAvoidLock)
{
    if(!canAvoidLock)pthread_mutex_lock(&queueLock);
#ifdef PERILLA_DEBUG
    if(n==perilla::MSG_QUEUE_DEFAULT_MAXSIZE){
        printf("Failed to Enqueue: Queue Overflow\n");
        exit(0);
    }
#endif
    buffer[rear] = package;
    prear = rear;
    rear = (rear+1)%perilla::MSG_QUEUE_DEFAULT_MAXSIZE;
    n++;
    if(!canAvoidLock)pthread_mutex_unlock(&queueLock);
}

Package* PackageQueue::dequeue(void)
{
    Package* package = 0;
    pthread_mutex_lock(&queueLock);
#ifdef PERILLA_DEBUG
    if(n<0){
        printf("Failed to Dequeue: Queue Empty\n");
        exit(0);
    }
#endif
    package = buffer[front];
    front = (front+1)%perilla::MSG_QUEUE_DEFAULT_MAXSIZE;
    n--;
    pthread_mutex_unlock(&queueLock);
    return package;
}

Package* PackageQueue::dequeue(bool canAvoidLock)
{
    Package* package = 0;
    if(!canAvoidLock)pthread_mutex_lock(&queueLock);
#ifdef PERILLA_DEBUG
    if(n<0){
        printf("Failed to Dequeue: Queue Empty\n");
        exit(0);
    }
#endif
    package = buffer[front];
    front = (front+1)%perilla::MSG_QUEUE_DEFAULT_MAXSIZE;
    n--;
    if(!canAvoidLock)pthread_mutex_unlock(&queueLock);
    return package;
}

Package* PackageQueue::getRear(void)
{
    Package* package = 0;
    pthread_mutex_lock(&queueLock);
    if(n) package = buffer[prear];
    pthread_mutex_unlock(&queueLock);
    return package;
}

Package* PackageQueue::getRear(bool canAvoidLock)
{
    Package* package = 0;
    if(!canAvoidLock)pthread_mutex_lock(&queueLock);
    if(n) package = buffer[prear];
    if(!canAvoidLock)pthread_mutex_unlock(&queueLock);
    return package;
}

Package* PackageQueue::getFront(void)
{
    Package* package = 0;
    pthread_mutex_lock(&queueLock);
    if(n) package = buffer[front];
    pthread_mutex_unlock(&queueLock);
    return package;
}

Package* PackageQueue::getFront(bool canAvoidLock)
{
    Package* package = 0;
    if(!canAvoidLock) pthread_mutex_lock(&queueLock);
    if(n) package = buffer[front];
    if(!canAvoidLock) pthread_mutex_unlock(&queueLock);
    return package;
}

void PackageQueue::emptyQueue(bool canAvoidLock){
    if(!canAvoidLock) pthread_mutex_lock(&queueLock);
    while(n){
	Package* p= dequeue(true);
	delete p;
    }
    if(!canAvoidLock) pthread_mutex_unlock(&queueLock);
}

PackageQueue::~PackageQueue()
{
    emptyQueue(true);   
}
