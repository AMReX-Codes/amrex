#ifndef _BL_WORKQUEUE_H_
#define _BL_WORKQUEUE_H_

//
// $Id: WorkQueue.H,v 1.3 2001-07-23 19:34:23 lijewski Exp $
//

#include <Thread.H>
#include <queue>

extern "C" void* WorkQueue_server(void*);
class WorkQueue;
class TimedWorkQueue;

class WorkQueue
{
public:
    friend void* WorkQueue_server(void*);
    class task;
    explicit WorkQueue(int maxthreads_ = 1 );
    virtual ~WorkQueue();
    int max_threads(int maxthreads_);
    int max_threads() const;
    int num_threads() const;
    void add(task* data);
    void drain();
    void wait();
    virtual void* server();
protected:
    ConditionVariable cv;
    Gate gate;
    std::queue<task*> wrkq;
    bool quit;
    bool eof;
    int maxthreads;
    int numthreads;
    int idlethreads;
    int tasks;
    void* operator new(size_t size);
    WorkQueue(const WorkQueue&);
    WorkQueue& operator=(const WorkQueue&);
};

class TimedWorkQueue
    : public WorkQueue
{
public:
    explicit TimedWorkQueue(int threads = 1, double timeout_ = -1.0 );
    double timeOut(double timeout);
    double timeOut() const;
private:
    double timeout;
    virtual void* server();
};

class WorkQueue::task
{
public:
    virtual void run() = 0;
    virtual ~task();
};

#endif /*_BL_WORKQUEUE_H_*/
