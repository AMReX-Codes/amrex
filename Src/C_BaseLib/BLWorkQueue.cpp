//
// $Id: BLWorkQueue.cpp,v 1.6 2001-09-14 21:03:06 lijewski Exp $
//

#include <winstd.H>
#include <cstdlib>
#include <cstdio>
#include <memory>
#include <queue>

#include <ParmParse.H>
#include <Thread.H>
#include <Profiler.H>
#include <WorkQueue.H>

namespace
{
Mutex print_mutex;
WorkQueue* bl_wrkq = 0;
int verbose = 0;
}

namespace BoxLib
{

WorkQueue&
theWorkQueue()
{
    return *bl_wrkq;
}
}

void
WorkQueue::Initialize ()
{
    ParmParse pp("workqueue");

    int maxthreads = 0;

    pp.query("maxthreads", maxthreads);
    pp.query("verbose", verbose);

    std::cout << "workqueue.maxthreads = " << maxthreads << std::endl;

    bl_wrkq = new WorkQueue(maxthreads);
}

void
WorkQueue::Finalize ()
{
    delete bl_wrkq;
}

#ifdef BL_THREADS
#define DPRINTF(arg)							\
do									\
{									\
    if ( verbose > 2 )							\
    {									\
	Lock<Mutex> lock(print_mutex);					\
	std::cout << "tid(" << Thread::getID() << "): "		\
		  << arg << std::endl;					\
    }									\
}									\
while (false)
#else
#define DPRINTF(arg)							\
do									\
{									\
    if ( verbose > 2 )							\
    {									\
	Lock<Mutex> lock(print_mutex);					\
	std::cout << "tid(" << 0 << "): "				\
		  << arg << std::endl;					\
    }									\
}									\
while (false)
#endif

WorkQueue::WorkQueue(int maxthreads_)
    : quit(false), eof(false), maxthreads(maxthreads_), numthreads(0), idlethreads(0), tasks(0)
{
    if ( maxthreads_ >= Thread::max_threads() )
    {
	BoxLib::Error("maxthreads_ in workqueue exceeds system limit");
    }
    if ( maxthreads_ < 0 )
    {
	BoxLib::Error("maxthreads_ must be >= 0");
    }
}

int
WorkQueue::max_threads() const
{
    return maxthreads;
}

int
WorkQueue::num_threads() const
{
    return numthreads;
}

WorkQueue::task::~task()
{
}

extern "C"
void*
WorkQueue_server(void* arg)
{
    BL_PROFILE("WorkQueue_server()");
    WorkQueue* wq = static_cast<WorkQueue*>(arg);
    return wq->server();
}

void*
WorkQueue::server()
{
    DPRINTF("A worker is starting");
    Lock<ConditionVariable> lock(cv);
    DPRINTF("Worker locked 0");
    for (;;)
    {
	if ( tasks == 0 && eof )
	{
	    gate.open(); gate.release();
	    eof = false;
	}
	DPRINTF("Worker waiting for work");
	while ( wrkq.empty() && !quit )
	{
	    idlethreads++;
	    cv.wait();
	    idlethreads--;
	}
	DPRINTF("Work queue: wrkq.empty()(" << wrkq.empty()<< "), "
		<< "quit(" << quit << "), "
		<< "eof("  << eof << "), "
		<< "tasks(" << tasks << ")");
	if ( !wrkq.empty() )
	{
	    std::auto_ptr<task> we(wrkq.front());
	    wrkq.pop();
	    if ( we.get() )
	    {
		eof = false;
		DPRINTF("Worker calling engine");
		cv.unlock();
		we->run();
		cv.lock();
		DPRINTF("Worker returning engine");
	    }
	    else
	    {
		DPRINTF("EOF reached");
		eof = true;
	    }
	    tasks--;
        }
	if ( wrkq.empty() && quit )
	{
	    DPRINTF("Worker shutting down");
	    if ( --numthreads == 0 )
	    {
		cv.broadcast();	// FIXME same predicate!
	    }
	    break;
        }
    }
    DPRINTF("Worker exiting");
    return 0;
}

void
WorkQueue::drain()
{
    Lock<ConditionVariable> lock(cv);
    if ( numthreads > 0 )
    {
	quit = true;
	if ( idlethreads > 0 )
	{
	    cv.broadcast();
        }
	while ( numthreads > 0 )
	{
	    cv.wait();		// FIXME??? using same cv for two predicates
        }
    }
}

WorkQueue::~WorkQueue()
{
    drain();
}

void
WorkQueue::add(task* item)
{
    std::auto_ptr<task> tsk(item);
    Lock<ConditionVariable> lock(cv);
    if ( maxthreads == 0 )
    {
	DPRINTF("maxthreads ==0 task");
	if ( item )
	{
	    item->run();
	}
	return;
    }
    wrkq.push(tsk.release());
    tasks++;
    if ( idlethreads > 0 )
    {
	DPRINTF("Signaling idle worker");
	cv.signal();
    }
    else if ( numthreads < maxthreads )
    {
	DPRINTF("Creating new worker");
	FunctionThread ft(WorkQueue_server, this, FunctionThread::Detached);
	numthreads++;
    }
}

void
WorkQueue::wait()
{
    add( 0 );
    if ( maxthreads )
    {
	gate.wait();
    }
    gate.close();
    DPRINTF("wait: finished...");
}

