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
    omp_init_lock(&packageLock);
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
    omp_init_lock(&packageLock);
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
    omp_init_lock(&packageLock);
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
    omp_set_lock(&packageLock);
    completed = true;
    omp_unset_lock(&packageLock);
}

void Package::completeRequest(bool lockIgnore)
{
    if(!lockIgnore)omp_set_lock(&packageLock);
    completed = true;
    if(!lockIgnore)omp_unset_lock(&packageLock);
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
    omp_init_lock(&packageLock);
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
    omp_init_lock(&queueLock);
}

int PackageQueue::queueSize(void)
{
    int size;
    omp_set_lock(&queueLock);
    size = n;
    omp_unset_lock(&queueLock);
    return size;
}

int PackageQueue::queueSize(bool lockIgnore)
{
    int size;
    if(!lockIgnore)omp_set_lock(&queueLock);
    size = n;
    if(!lockIgnore)omp_unset_lock(&queueLock);
    return size;
}

void PackageQueue::enqueue(Package* package)
{
    omp_set_lock(&queueLock);
    buffer[rear] = package;
    prear = rear;    
    rear = (rear+1)%perilla::MSG_QUEUE_DEFAULT_MAXSIZE;
    n++;
    omp_unset_lock(&queueLock);
}

void PackageQueue::enqueue(Package* package, bool lockIgnore)
{
    if(!lockIgnore)omp_set_lock(&queueLock);
    buffer[rear] = package;
    prear = rear;
    rear = (rear+1)%perilla::MSG_QUEUE_DEFAULT_MAXSIZE;
    n++;
    if(!lockIgnore)omp_unset_lock(&queueLock);
}

Package* PackageQueue::dequeue(void)
{
    Package* package = 0;
    omp_set_lock(&queueLock);
    package = buffer[front];
    front = (front+1)%perilla::MSG_QUEUE_DEFAULT_MAXSIZE;
    n--;
    omp_unset_lock(&queueLock);
    return package;
}

Package* PackageQueue::dequeue(bool lockIgnore)
{
    lockIgnore = false;
    Package* package = 0;
    if(!lockIgnore)omp_set_lock(&queueLock);
    if(n<0)
	std::cout<< "Q size " << n << " front " << front <<std::endl;
    package = buffer[front];
    front = (front+1)%perilla::MSG_QUEUE_DEFAULT_MAXSIZE;
    n--;
    if(!lockIgnore)omp_unset_lock(&queueLock);
    return package;
}

Package* PackageQueue::getRear(void)
{
    Package* package = 0;
    omp_set_lock(&queueLock);
    package = buffer[prear];
    omp_unset_lock(&queueLock);
    return package;
}
Package* PackageQueue::getRear(bool lockIgnore)
{
    Package* package = 0;
    if(!lockIgnore)omp_set_lock(&queueLock);
    package = buffer[prear];
    if(!lockIgnore)omp_unset_lock(&queueLock);
    return package;
}

Package* PackageQueue::getFront(void)
{
    Package* package = 0;
    omp_set_lock(&queueLock);
    package = buffer[front];
    omp_unset_lock(&queueLock);
    return package;
}

Package* PackageQueue::getFront(bool lockIgnore)
{
    Package* package = 0;
    if(!lockIgnore) omp_set_lock(&queueLock);
    package = buffer[front];
    if(!lockIgnore) omp_unset_lock(&queueLock);
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

