#include <Thread.H>

#include <unistd.h>
#include <sys/resource.h>
#include <sys/time.h>

#include <stdexcept>
#include <string>
#include <iostream>
#include <limits>

#include <cstdio>
#include <ctime>
#include <cerrno>
#include <cstdlib>

#define THREAD_ASSERT(x)						\
do									\
{									\
  if ( !(x) )								\
    {									\
      BoxLib::Error(#x );	\
    }									\
}									\
while ( false )

int Thread::m_next_id = 0;
Mutex Thread::next_id;
ThreadSpecificData<Thread> Thread::m_self(0,0);

namespace
{
std::string
the_message_string(const char* file, int line, const char* call, int status)
{
    //
    // Should be large enough.
    //
    const int DIM = 1024;
    char buf[DIM];
    if ( status )
    {
	std::sprintf(buf, "File %s, line %d, %s: %s", file, line, call, std::strerror(status));
    }
    else
    {
	std::sprintf(buf, "File %s, line %d, %s", file, line, call);
    }
    buf[DIM-1] = '\0';		// Just to be safe.
    return std::string(buf);
}
}

Thread::thread_error::thread_error(const char* file, int line, const char* call, int status)
    : std::runtime_error( the_message_string(file, line, call, status) )
{
}

int
Thread::get_next_id()
{
    Lock<Mutex> lock(next_id);
    int nid = m_next_id++;
    return nid;
}

Thread::Thread()
    : m_status(Created), m_id(get_next_id()), m_jod(false)
{
}

pthread_t&
Thread::theThread()
{
    return m_tid;
}

int
Thread::getID() const
{
    return m_id;
}

void*
Thread::_doit(void* arg)
{
    try
    {
	m_self.set(static_cast<Thread*>(arg));
	return static_cast<Thread*>(arg)->work();
    }
    catch( Thread::thread_error& bad )
    {
	std::cerr << "_doit caught " << bad.what() << std::endl << std::flush;
	std::abort();
    }
    catch(...)
    {
	// catch all uncaught exceptions from threads
	std::cerr << "Fatal...uncaught exception in _doit" << std::endl;
	std::abort();
    }
    return 0;
}

Thread*
Thread::self()
{
    return m_self.get();
}

void
Thread::sleep(const Time& spec_)
{
#ifndef BL_USE_NANOSLEEP
    unsigned long tosleep = spec_.as_long();
    if ( tosleep > std::numeric_limits<unsigned long>::max()/1000000 )
    {
	::sleep(tosleep);
    }
    else
    {
	::usleep(tosleep*1000000);
    }
#else
    Time spec = spec_;
    Time rem;
    while ( int status = ::nanosleep(&spec, &rem) )
    {
	if ( errno != EINVAL )
	{
	    throw thread_error(__FILE__, __LINE__, "nanosleep", errno);
	}
	spec = rem;
    }
#endif
}


//
// Barrier
//

Barrier::Barrier(int i)
    : count(0), n_sleepers(0), releasing(false)
{
    init(i);
}

void
Barrier::init(int i)
{
    THREAD_ASSERT( !releasing );
    THREAD_ASSERT( n_sleepers == 0 );
    count = i;
}

void
Barrier::wait()
{
    BL_PROFILE( BL_PROFILE_THIS_NAME() + "::wait()" );
    bool release = false;
    lock();
    // If previous cycle still releasing, wait
    // THREAD_ASSERT ( !releasing );
    while ( releasing )
    {
	ConditionVariable::wait();
    }
    if ( ++n_sleepers == count )
    {
	release = releasing = true;
    }
    else
    {
	// A poor thread cancelation Site
	Thread::CancelState tmp = Thread::setCancelState(Thread::Disable);
	while ( !releasing )
	{
	    ConditionVariable::wait();
	}
	Thread::setCancelState(tmp);
    }
    if ( --n_sleepers == 0 )
    {
	releasing = false;
	release = true;             // Wake up waiters (if any) for next cycle
    }
    unlock();
    if ( release )
    {
	broadcast();
    }
}


//
// Semaphore
//

Semaphore::Semaphore(int val_)
    : value(val_)
{
}

void
Semaphore::wait()
{
    BL_PROFILE( BL_PROFILE_THIS_NAME() + "::wait()" );
    lock();
    while ( value == 0 )
    {
	ConditionVariable::wait();
    }
    value--;
    unlock();
}

bool
Semaphore::trywait()
{
    lock();
    if ( value == 0 )
    {
	unlock();
	return false;
    }
    value--;
    unlock();
    return true;
}

void
Semaphore::post()
{
    lock();
    value++;
    unlock();
    signal();
}


//
//
//

SemaphoreB::SemaphoreB(int val_)
    : val(val_)
{}

int
SemaphoreB::down()
{
    lock();
    while (val <= 0)
    {
	wait();
    }
    int t = --val;
    unlock();

    return t;
}

int
SemaphoreB::up()
{
    lock();
    int t = ++val;
    unlock();
    signal();
    return t;
}

int
SemaphoreB::decrement()
{
    lock();
    int t = --val;
    unlock();
    return t;
}

int
SemaphoreB::value()
{
    lock();
    int t = val;
    unlock();
    return t;
}

//
// SingleBarrier
//

SingleBarrier::SingleBarrier(int i)
    : count(i), n_posters(0), n_waiters(0), releasing(false)
{
}

void
SingleBarrier::wait()
{
    bool release = false;
    lock();
    n_waiters++;
    while ( !releasing )
    {
	ConditionVariable::wait();
    }
    if ( --n_waiters == 0 )
    {
	releasing = false;
	release = true;             // Wake up waiters (if any) for next cycle
	n_posters=0;
    }
    unlock();
    if ( release )
    {
	broadcast();
    }
}

void
SingleBarrier::post()
{
    bool release = false;
    lock();
    // If previous cycle still releasing, wait
    while ( releasing )
    {
	ConditionVariable::wait();
    }
    if ( ++n_posters == count )
    {
	releasing = true;
	release = true;             // Wake up waiters (if any) for next cycle
    }
    unlock();
    if ( release )
    {
	broadcast();
    }
}


//
//Gate
//

Gate::Gate()
    : closed(true)
{
}

void
Gate::open()
{
    lock();
    closed = false;
    broadcast();
    unlock();
}

void
Gate::close()
{
    lock();
    closed = true;
    unlock();
}

void
Gate::release()
{
    broadcast();
}

void
Gate::wait()
{
    BL_PROFILE( BL_PROFILE_THIS_NAME() + "::wait()" );
    lock();
    while ( closed )
    {
	ConditionVariable::wait();
    }
    unlock();
}


//
// Lock<Semaphore> specialization
//

Lock<Semaphore>::Lock(Semaphore& sem_)
    : sem(sem_)
{
    sem.wait();
}

Lock<Semaphore>::~Lock()
{
    sem.post();
}

#ifdef BL_PTHREADS

namespace
{
extern "C"
{
    typedef void* (*thr_vpvp)(void*);
    typedef void (*thr_vvp)(void*);
}
}



//
// Thread
//

Thread::~Thread()
{
    if  ( m_status == Running && !m_jod )
    {
	THREAD_REQUIRE( pthread_detach(m_tid) );
    }
}

int
Thread::max_threads()
{
#ifdef PTHREAD_THREADS_MAX
    return PTHREAD_THREADS_MAX;
#else
    return 64;
#endif
}

void
Thread::run(DetachState instate)
{
    run(Attr(instate));
}

void
Thread::run(const Attr& attr)
{
    THREAD_ASSERT( m_status == Created );
    THREAD_ASSERT( !m_jod );
    if ( !attr.isJoinable() )
    {
	m_jod = true;
    }
    THREAD_REQUIRE( pthread_create(&m_tid, &(attr.m_attr), reinterpret_cast<thr_vpvp>(_doit), this) );
    m_status = Running;
}

void
Thread::detach()
{
    THREAD_ASSERT( m_status == Running );
    THREAD_ASSERT( !m_jod );
    THREAD_REQUIRE( pthread_detach(m_tid) );
    m_jod =true;
}

void
Thread::exit(void* st)
{
    if ( Thread* me = self() )
    {
	me->m_mutex.lock();
	me->m_status = Stopped;
	me->m_mutex.unlock();
	// if ( me->m_jod ) delete me;
    }
    pthread_exit(st);
}

void*
Thread::join() const
{
    THREAD_ASSERT( m_status == Running );
    THREAD_ASSERT( !m_jod );
    void* exitst;
    THREAD_REQUIRE( pthread_join(m_tid, &exitst) );
    m_jod = true;
    return exitst;
}

void
Thread::yield()
{
#ifdef _POSIX_PRIORITY_SCHEDULING
    sched_yield();
#endif
}

Thread::CancelState
Thread::setCancelState(CancelState state)
{
    CancelState result;
    int newstate;
    switch ( state )
    {
    case Enable:
	newstate = PTHREAD_CANCEL_ENABLE;
	break;
    case Disable:
	newstate = PTHREAD_CANCEL_DISABLE;
	break;
    }
    int oldstate;
    THREAD_REQUIRE( pthread_setcancelstate(newstate, &oldstate) );
    switch ( oldstate )
    {
    case PTHREAD_CANCEL_ENABLE:
	result = Enable;
	break;
    case PTHREAD_CANCEL_DISABLE:
	result = Disable;
	break;
    }
    return result;
}

Thread::CancelType
Thread::setCancelType(CancelType t)
{
    CancelType result;
    int newtype;
    switch ( t )
    {
    case Deferred:
	newtype = PTHREAD_CANCEL_DEFERRED;
	break;
    case Asynchronous:
	newtype = PTHREAD_CANCEL_ASYNCHRONOUS;
	break;
    }
    int oldtype;
    THREAD_REQUIRE( pthread_setcanceltype(newtype, &oldtype) );
    switch ( oldtype )
    {
    case PTHREAD_CANCEL_DEFERRED:
	result = Deferred;
	break;
    case PTHREAD_CANCEL_ASYNCHRONOUS:
	result = Asynchronous;
	break;
    }
    return result;
}

void
Thread::testcancel()
{
    pthread_testcancel();
}

void
Thread::cancel()
{
    THREAD_ASSERT( m_status == Running );
    THREAD_REQUIRE( pthread_cancel(m_tid) );
}


//
// Thread::Attr s.
//

bool
Thread::Attr::isJoinable() const
{
    return m_is_joinable;
}

Thread::Attr::Attr(DetachState inistate)
    : m_is_joinable( true )
{
    THREAD_REQUIRE( pthread_attr_init(&m_attr) );
    if ( Debug::Ndebug )
    {
	int l;
	THREAD_REQUIRE( pthread_attr_getdetachstate(&m_attr, &l) );
	THREAD_ASSERT(l == PTHREAD_CREATE_JOINABLE);
    }
    if ( inistate == Detached )
    {
	m_is_joinable = false;
	set_detachstate(inistate);
    }
}

Thread::Attr::~Attr()
{
    THREAD_REQUIRE( pthread_attr_destroy(&m_attr) );
}

Thread::DetachState
Thread::Attr::get_detachstate() const
{
    int l;
    THREAD_REQUIRE( pthread_attr_getdetachstate(&m_attr, &l) );
    DetachState result;
    switch ( l )
    {
    case PTHREAD_CREATE_JOINABLE:
	result = Joinable;
	break;
    case PTHREAD_CREATE_DETACHED:
	result = Detached;
	break;
    }
    return result;
}

Thread::DetachState
Thread::Attr::set_detachstate(Thread::DetachState ll)
{
    int detachstate;
    switch ( ll )
    {
    case Joinable:
	detachstate = PTHREAD_CREATE_JOINABLE;
	break;
    case Detached:
	detachstate = PTHREAD_CREATE_DETACHED;
	break;
    }
    Thread::DetachState result = get_detachstate();
    THREAD_REQUIRE( pthread_attr_setdetachstate(&m_attr, detachstate) );
    return result;
}

#ifdef _POSIX_THREAD_PRIORITY_SCHEDULING
Thread::Attr::SchedPolicy
Thread::Attr::get_schedpolicy() const
{
    int l;
    THREAD_REQUIRE( pthread_attr_getschedpolicy(&m_attr, &l) );
    SchedPolicy result;
    switch ( l )
    {
    case SCHED_OTHER:
	result = Regular;
	break;
    case SCHED_RR:
	result = RoundRobin;
	break;
    case SCHED_FIFO:
	result = FiFo;
	break;
    }
    return result;
}

Thread::Attr::SchedPolicy
Thread::Attr::set_schedpolicy(SchedPolicy schedpolicy)
{
    int ll;
    switch ( schedpolicy )
    {
    case Regular: ll = SCHED_OTHER; break;
    case RoundRobin: ll = SCHED_RR; break;
    case FiFo: ll = SCHED_FIFO; break;
    }
    SchedPolicy l = get_schedpolicy();
    THREAD_REQUIRE( pthread_attr_setschedpolicy(&m_attr, ll) );
    return l;
}

Thread::Attr::InheritSched
Thread::Attr::get_inheritsched() const
{
    int l;
    THREAD_REQUIRE( pthread_attr_getinheritsched(&m_attr, &l) );
    InheritSched result;
    switch ( l )
    {
    case PTHREAD_EXPLICIT_SCHED:
	result = Explicit;
	break;
    case PTHREAD_INHERIT_SCHED:
	result = Inherit;
	break;
    }
    return result;
}

Thread::Attr::InheritSched
Thread::Attr::set_inheritsched(InheritSched inheritsched)
{
    int ll;
    switch ( inheritsched )
    {
    case Explicit: ll = PTHREAD_EXPLICIT_SCHED; break;
    case Inherit: ll = PTHREAD_INHERIT_SCHED; break;
    }
    InheritSched l = get_inheritsched();
    THREAD_REQUIRE( pthread_attr_setinheritsched(&m_attr, ll) );
    return l;
}

Thread::Attr::SchedScope
Thread::Attr::get_scope() const
{
    int l;
    THREAD_REQUIRE( pthread_attr_getscope(&m_attr, &l) );
    SchedScope result;
    switch ( l )
    {
    case PTHREAD_SCOPE_SYSTEM:
	result = System;
	break;
    case PTHREAD_SCOPE_PROCESS:
	result = Process;
	break;
    }
    return result;
}

Thread::Attr::SchedScope
Thread::Attr::set_scope(SchedScope scope)
{
    int ll;
    switch ( scope )
    {
    case System: ll = PTHREAD_SCOPE_SYSTEM; break;
    case Process: ll = PTHREAD_SCOPE_PROCESS; break;
    }
    SchedScope l = get_scope();
    THREAD_REQUIRE( pthread_attr_setscope(&m_attr, ll) );
    return l;
}

sched_param
Thread::Attr::get_schedparam() const
{
    sched_param l;
    THREAD_REQUIRE( pthread_attr_getschedparam(&m_attr, &l) );
    return l;
}

sched_param
Thread::Attr::set_schedparam(const sched_param& schedparam)
{
    sched_param l = get_schedparam();
    THREAD_REQUIRE( pthread_attr_setschedparam(&m_attr, &schedparam) );
    return l;
}

#endif


//
//Mutex
//

Mutex::Mutex()
{
    THREAD_REQUIRE( pthread_mutex_init(&m_mutex, 0) );
}

Mutex::Mutex(const Attr& attr)
{
#ifdef _AIX_PTHREADS_D7
    THREAD_REQUIRE( pthread_mutex_init(&m_mutex,
				       const_cast<pthread_mutexattr_t*>(&attr.attribute())) );
#else
    THREAD_REQUIRE( pthread_mutex_init(&m_mutex, &attr.attribute()) );
#endif
}

Mutex::~Mutex()
{
    THREAD_REQUIRE( pthread_mutex_destroy(&m_mutex) );
}

void
Mutex::lock()
{
    THREAD_REQUIRE( pthread_mutex_lock(&m_mutex) );
}

pthread_mutex_t&
Mutex::theMutex()
{
    return m_mutex;
}

bool
Mutex::trylock()
{
    int status = pthread_mutex_trylock(&m_mutex);
    if ( status == 0 ) return true;
    if ( status == EBUSY ) return false;
    throw Thread::thread_error(__FILE__,__LINE__,"pthread_mutex_trylock(&m_mutex", status);
}

void
Mutex::unlock()
{
    THREAD_REQUIRE( pthread_mutex_unlock(&m_mutex) );
}


//
// Mutex Atrributes
//

Mutex::Attr::Attr()
{
    THREAD_REQUIRE( pthread_mutexattr_init(&m_attr) );
}

Mutex::Attr::~Attr()
{
    THREAD_REQUIRE( pthread_mutexattr_destroy(&m_attr) );
}

const pthread_mutexattr_t&
Mutex::Attr::attribute() const
{
    return m_attr;
}

pthread_mutexattr_t&
Mutex::Attr::attribute()
{
    return m_attr;
}


//
// ConditionVariable
//

ConditionVariable::ConditionVariable()
{
    THREAD_REQUIRE( pthread_cond_init(&m_cv, 0) );
}

ConditionVariable::ConditionVariable(const Attr& attr)
{
#ifdef _AIX_PTHREADS_D7
    THREAD_REQUIRE( pthread_cond_init(&m_cv,
				      const_cast<pthread_condattr_t*>(&attr.attribute())) );
#else
    THREAD_REQUIRE( pthread_cond_init(&m_cv, &attr.attribute()) );
#endif
}

ConditionVariable::~ConditionVariable()
{
    THREAD_REQUIRE( pthread_cond_destroy(&m_cv) );
}

void
ConditionVariable::signal()
{
    THREAD_REQUIRE( pthread_cond_signal(&m_cv) );
}

void
ConditionVariable::broadcast()
{
    THREAD_REQUIRE( pthread_cond_broadcast(&m_cv) );
}

void
ConditionVariable::wait()
{
    THREAD_REQUIRE( pthread_cond_wait(&m_cv, &theMutex()) );
}

void
ConditionVariable::wait(Mutex& m)
{
    THREAD_REQUIRE( pthread_cond_wait(&m_cv, &m.theMutex()) );
}

bool
ConditionVariable::timedwait(const Time& abstime)
{
    int status = pthread_cond_timedwait(&m_cv, &theMutex(), &abstime);

    if ( status == 0 ) return false;
    if ( status == ETIMEDOUT ) return true;
    throw Thread::thread_error(__FILE__, __LINE__, "pthread_cond_timedwait(&m_cv, &theMutex(), &abstime)", status);
}

bool
ConditionVariable::timedwait(Mutex& m, const Time& abstime)
{
    int status = pthread_cond_timedwait(&m_cv, &m.theMutex(), &abstime);

    if ( status == 0 ) return false;
    if ( status == ETIMEDOUT ) return true;
    throw Thread::thread_error(__FILE__, __LINE__, "pthread_cond_timedwait", status);
}

pthread_cond_t&
ConditionVariable::theCV()
{
    return m_cv;
}

//
// ConditionVariable Attributes
//

ConditionVariable::Attr::Attr()
{
    THREAD_REQUIRE( pthread_condattr_init(&m_attr) );
}

ConditionVariable::Attr::~Attr()
{
    THREAD_REQUIRE( pthread_condattr_destroy(&m_attr) );
}

const pthread_condattr_t&
ConditionVariable::Attr::attribute() const
{
    return m_attr;
}

pthread_condattr_t&
ConditionVariable::Attr::attribute()
{
    return m_attr;
}

#ifdef _POSIX_THREAD_PROCESS_SHARED
ConditionVariable::Attr::ProcessShared
ConditionVariable::Attr::get_pshared() const
{
    ProcessShared result;
    int l;
    THREAD_REQUIRE( pthread_condattr_getpshared(&m_attr, &l) );
    switch ( l )
    {
    case PTHREAD_PROCESS_SHARED:
	result = Shared;
	break;
    case PTHREAD_PROCESS_PRIVATE:
	result = Private;
	break;
    }
    return result;
}

ConditionVariable::Attr::ProcessShared
ConditionVariable::Attr::set_pshared(ProcessShared v)
{
    int ll;
    switch ( v )
    {
    case Shared:
	ll = PTHREAD_PROCESS_SHARED;
	break;
    case Private:
	ll = PTHREAD_PROCESS_PRIVATE;
	break;
    }
    ProcessShared result = get_pshared();
    THREAD_REQUIRE( pthread_condattr_setpshared(&m_attr, ll) );
    return result;
}
#endif


//
// RecursiveMutex
//

RecursiveMutex::RecursiveMutex()
    : count(0), owned(false)
{
}

RecursiveMutex::~RecursiveMutex()
{
}

void
RecursiveMutex::lock()
{
    const pthread_t tid = pthread_self();
    Mutex::lock();
    while ( owned && !pthread_equal(m_tid, tid) )
    {
	wait();
    }
    m_tid = tid;
    owned = true;
    count++;
    Mutex::unlock();
}

bool
RecursiveMutex::trylock()
{
    const pthread_t tid = pthread_self();
    Mutex::lock();
    if ( !owned || pthread_equal(m_tid, tid) )
    {
	m_tid = tid;
	owned = true;
	count++;
	Mutex::unlock();
	return true;
    }
    Mutex::unlock();
    return false;
}

void
RecursiveMutex::unlock()
{
    Mutex::lock();
    if ( --count == 0 )
    {
	owned = false;
	Mutex::unlock();
	signal();
	return;
    }
    Mutex::unlock();
}


//
//Thread Specific Data
//

ThreadSpecificData<void>::ThreadSpecificData(void (*tsd)(void*))
{
    THREAD_REQUIRE( pthread_key_create(&key, reinterpret_cast<thr_vvp>(tsd)) );
    THREAD_ASSERT(get() == 0);
}

ThreadSpecificData<void>::~ThreadSpecificData()
{
    THREAD_REQUIRE( pthread_key_delete(key) );
}

void*
ThreadSpecificData<void>::set(const void* v)
{
    void* ov = pthread_getspecific(key);
    THREAD_REQUIRE( pthread_setspecific(key, v) );
    return ov;
}

void*
ThreadSpecificData<void>::get() const
{
    return pthread_getspecific(key);
}


//
// RWLock
//

void
RWLock::_rdlock_cleanup(void *arg)
{
    static_cast<ConditionVariable*>(arg)->unlock();
}

void
RWLock::_wrlock_cleanup(void *arg)
{
    RWLock* rwp = static_cast<RWLock*>(arg);
    // Was the only queued writer and lock is available for readers.
    // Called through cancellation clean-up so lock is held at entry.
    if ( (--rwp->waiters == 0) && (rwp->state >= 0) )
    {
	rwp->readers.broadcast();
    }
    rwp->readers.unlock();
}

RWLock::RWLock()
    : state(0), waiters(0)
{
}

RWLock::~RWLock()
{
}

void
RWLock::rdlock()
{
    readers.lock();
    pthread_cleanup_push(reinterpret_cast<thr_vvp>(_rdlock_cleanup), &readers);
    // active or queued writers
    while ( (state < 0) || waiters )
    {
	readers.wait();
    }
    state++;
    pthread_cleanup_pop(1);
}

bool
RWLock::tryrdlock()
{
    bool result = false;
    Lock<ConditionVariable> lock(readers);
    // available and no writers queued
    if ( (state >= 0) && !waiters )
    {
	state++;
	result = true;
    }
    return result;
}

void
RWLock::wrlock()
{
    readers.lock();
    waiters++;			// another writer queued
    pthread_cleanup_push(reinterpret_cast<thr_vvp>(_wrlock_cleanup), static_cast<void*>(this));

    while ( state )
    {
	writers.wait();
    }
    state = -1;
    pthread_cleanup_pop(1);
}

bool
RWLock::trywrlock()
{
    bool result = false;
    Lock<ConditionVariable> lock(readers);
    // no readers, no writers, no writers queued
    if ( !state && !waiters )
    {
	state = -1;
	result = true;
    }
    return result;
}

void
RWLock::unlock()
{
    Lock<ConditionVariable> lock(readers);
    if ( state == -1 )		// writer releasing
    {
	state = 0;		// mark as available
	if ( waiters )
	{
	    writers.signal();	// writers queued
	}
	else
	{
	    readers.broadcast();
	}
    }
    else
    {
	if ( --state == 0 )
	{
	    writers.signal();	// no more readers
	}
    }
}



FunctionThread::FunctionThread(Thread_Function func_, void* arg_, Thread::DetachState st)
    : m_jod(false)
{
    BL_PROFILE( BL_PROFILE_THIS_NAME() + "::FunctionThread()" );
    Thread::Attr attr(st);
    pthread_attr_t a = attr.theAttribute();
    if ( !attr.isJoinable() )
    {
	m_jod = true;
    }
    THREAD_REQUIRE( pthread_create(&m_tid, &a, func_, arg_) );
}

FunctionThread::~FunctionThread()
{
    detach();
}

void*
FunctionThread::join() const
{
    BL_PROFILE( BL_PROFILE_THIS_NAME() + "::join()" );
    void* ret;
    if ( !m_jod )
    {
	THREAD_REQUIRE( pthread_join(m_tid, &ret) );
	m_jod = true;
    }
    return ret;
}

void
FunctionThread::detach() const
{
    BL_PROFILE( BL_PROFILE_THIS_NAME() + "::detach()" );
    if ( !m_jod )
    {
	THREAD_REQUIRE( pthread_detach(m_tid) );
    }
    m_jod = true;
}

#else

void
Thread::run(DetachState inistate)
{
    _doit(this);
}

void
Thread::exit(void*)
{
    std::exit(0);
}

Thread::CancelState
Thread::setCancelState(CancelState)
{
    return Enable;
}


//
// A ConditionVariable
//

bool
ConditionVariable::timedwait(const Time& t)
{
    Thread::sleep(t);
    return true;
}

bool
ConditionVariable::timedwait(Mutex&, const Time& t)
{
    Thread::sleep(t);
    return true;
}


//
//
//
ThreadSpecificData<void>::ThreadSpecificData(void (*tsd_)(void*))
    : v(0), tsd(tsd_)
{
}

ThreadSpecificData<void>::~ThreadSpecificData()
{
    // std::cout << "got here" << std::endl;
    // (*tsd)(v);
}

void*
ThreadSpecificData<void>::set(const void* v_)
{
    return v = const_cast<void*>(v_);
}

void*
ThreadSpecificData<void>::get() const
{
    return v;
}

// FuctioinThread
FunctionThread::FunctionThread(Thread_Function func, void* arg_, Thread::DetachState st)
    : m_jod(false)
{
    func(arg_);
}

FunctionThread::~FunctionThread()
{
    detach();
}

void*
FunctionThread::join() const
{
    m_jod = true;
    return 0;
}

void
FunctionThread::detach() const
{
    m_jod = true;
}

#endif
