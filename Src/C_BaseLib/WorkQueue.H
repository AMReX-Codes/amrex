#ifndef _BL_WORKQUEUE_H_
#define _BL_WORKQUEUE_H_

//
// $Id: WorkQueue.H,v 1.5 2001-11-16 18:19:33 car Exp $
//

#include <Thread.H>
#include <queue>

extern "C" void* WorkQueue_server(void*);
class WorkQueue;

class WorkQueue
{
public:
    friend void* WorkQueue_server (void*);
    class task;
    explicit WorkQueue (int maxthreads_);
    virtual ~WorkQueue ();
    int max_threads () const;
    int num_threads () const;
    void add (task* data);
    void drain ();
    void wait ();
    virtual void* server ();
    static void Initialize ();
    static void Finalize ();
protected:
    mutable ConditionVariable cv;
    Gate gate;
    std::queue<task*> wrkq;
    bool quit;
    bool eof;
    const int maxthreads;
    int numthreads;
    int idlethreads;
    int tasks;
    WorkQueue (const WorkQueue&);
    WorkQueue& operator= (const WorkQueue&);
};

class WorkQueue::task
{
public:
    virtual void run () = 0;
    virtual ~task ();
};

namespace BoxLib
{
WorkQueue& theWorkQueue ();
}


#endif /*_BL_WORKQUEUE_H_*/
