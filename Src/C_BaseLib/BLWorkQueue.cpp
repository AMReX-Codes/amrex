//
// $Id: BLWorkQueue.cpp,v 1.3 2001-07-26 20:27:27 car Exp $
//

#include <Thread.H>

#include <cstdlib>
#include <cstdio>
#include <memory>
#include <queue>

#include <WorkQueue.H>

namespace
{
Mutex print_mutex;
}

//#define BL_WORKQUEUE_TRACE
#ifdef BL_WORKQUEUE_TRACE
#define DPRINTF(arg)							\
do									\
{									\
    Lock<Mutex> lock(print_mutex);					\
    std::cout << "tid(" << (long)(pthread_self()) << "): "		\
	      << arg << std::endl;					\
}									\
while (false)
#else
#define DPRINTF(arg)
#endif

int
WorkQueue::max_threads(int maxthreads_)
{
    int ot = maxthreads;
    maxthreads = maxthreads_;
    return ot;
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

WorkQueue::WorkQueue(int maxthreads_)
    : quit(false), eof(false), maxthreads(maxthreads_), numthreads(0), idlethreads(0), tasks(0)
{
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
	FunctionThread ft(WorkQueue_server, this, Thread::Detached);
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

TimedWorkQueue::TimedWorkQueue(int maxthreads_, double timeout_)
    : WorkQueue(maxthreads_), timeout(timeout_)
{
}

double
TimedWorkQueue::timeOut(double timeout_)
{
    double ot = timeout;
    timeout = timeout_;
    return ot;
}

double
TimedWorkQueue::timeOut() const
{
    return timeout;
}

void*
TimedWorkQueue::server()
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
	bool timedout = false;
	BoxLib::Time t = BoxLib::Time::get_time();
	t += timeout;
	while ( wrkq.empty() && !quit && !eof )
	{
	    idlethreads++;
	    bool status = cv.timedwait(t);
	    idlethreads--;
	    if ( status )
	    {
		DPRINTF("Worker wait timed out");
		timedout = true;
		break;
	    }
	}
	DPRINTF("Work queue: wrkq.empty()(" << wrkq.empty() << "), "
		<< "quit(" << quit << "), "
		<< "eof(" << eof << ")");
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
	if ( wrkq.empty() && timedout )
	{
	    DPRINTF("engine terminating due to timeout.");
	    numthreads--;
	    break;
        }
    }
    DPRINTF("Worker exiting");
    return 0;
}
