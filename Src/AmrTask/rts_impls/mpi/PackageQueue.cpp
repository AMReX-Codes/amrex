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
    notified = false;
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
    notified = false;
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
#ifdef PERILLA_DEBUG
    memcheck.add(memcheck.genKey(this), (void*)this, "Package");
#endif
}

Package::Package(int src, int dest, int size)
{
    source = src;
    destination = dest;
    databuf = new double[size];
    bufSize = size;
    source = 0;
    destination = 0;
    completed = false;
    notified = false;
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

void Package::completeRequest(bool lockIgnore)
{
    if(!lockIgnore)pthread_mutex_lock(&packageLock);
    completed = true;
    if(!lockIgnore)pthread_mutex_unlock(&packageLock);
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
    notified = false;
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

int PackageQueue::queueSize(bool lockIgnore)
{
    int size;
    if(!lockIgnore)pthread_mutex_lock(&queueLock);
    size = n;
    if(!lockIgnore)pthread_mutex_unlock(&queueLock);
    return size;
}

void PackageQueue::enqueue(Package* package)
{
    pthread_mutex_lock(&queueLock);
    buffer[rear] = package;
    prear = rear;    
    rear = (rear+1)%perilla::MSG_QUEUE_DEFAULT_MAXSIZE;
    n++;
    pthread_mutex_unlock(&queueLock);
}

void PackageQueue::enqueue(Package* package, bool lockIgnore)
{
    if(!lockIgnore)pthread_mutex_lock(&queueLock);
    buffer[rear] = package;
    prear = rear;
    rear = (rear+1)%perilla::MSG_QUEUE_DEFAULT_MAXSIZE;
    n++;
    if(!lockIgnore)pthread_mutex_unlock(&queueLock);
}

Package* PackageQueue::dequeue(void)
{
    Package* package = 0;
    pthread_mutex_lock(&queueLock);
    package = buffer[front];
    front = (front+1)%perilla::MSG_QUEUE_DEFAULT_MAXSIZE;
    n--;
    pthread_mutex_unlock(&queueLock);
    return package;
}

Package* PackageQueue::dequeue(bool lockIgnore)
{
    lockIgnore = false;
    Package* package = 0;
    if(!lockIgnore)pthread_mutex_lock(&queueLock);
    if(n<0)
	std::cout<< "Q size " << n << " front " << front <<std::endl;
    package = buffer[front];
    front = (front+1)%perilla::MSG_QUEUE_DEFAULT_MAXSIZE;
    n--;
    if(!lockIgnore)pthread_mutex_unlock(&queueLock);
    return package;
}

Package* PackageQueue::getRear(void)
{
    Package* package = 0;
    pthread_mutex_lock(&queueLock);
    package = buffer[prear];
    pthread_mutex_unlock(&queueLock);
    return package;
}
Package* PackageQueue::getRear(bool lockIgnore)
{
    Package* package = 0;
    if(!lockIgnore)pthread_mutex_lock(&queueLock);
    package = buffer[prear];
    if(!lockIgnore)pthread_mutex_unlock(&queueLock);
    return package;
}

Package* PackageQueue::getFront(void)
{
    Package* package = 0;
    pthread_mutex_lock(&queueLock);
    package = buffer[front];
    pthread_mutex_unlock(&queueLock);
    return package;
}

Package* PackageQueue::getFront(bool lockIgnore)
{
    Package* package = 0;
    if(!lockIgnore) pthread_mutex_lock(&queueLock);
    package = buffer[front];
    if(!lockIgnore) pthread_mutex_unlock(&queueLock);
    return package;
}

void PackageQueue::emptyQueue(){
    while(n){
	Package* p= dequeue(false);
	delete p;
    }
}

PackageQueue::~PackageQueue()
{
    emptyQueue();   
}

