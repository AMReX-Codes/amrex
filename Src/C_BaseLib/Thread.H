#ifndef _BL_THREAD_H_
#define _BL_THREAD_H_

#ifdef BL_PTHREADS
#include <pthread.h>
#include <unistd.h>
#else
#if !( defined(__GNUC__) || defined(BL_AIX) )
typedef int pthread_mutex_t;
typedef int pthread_cond_t;
typedef int pthread_key_t;
typedef int pthread_t;
#endif
#endif

#include <stdexcept>
#include <Profiler.H>
#include <Utility.H>

//
// Generic exception class for POSIX thread errors.
//

//
// An abstract base class for defining threaded operations.
//
// You must be very careful when using Pthread objects to ensure
// that the Pthread object exists at least as long as the POSIX
// thread executing the work() routine is running.
// Most of this was written by Mike Lijewski <mjlijewski@lbl.gov>

extern "C"
{
    typedef void* (*Thread_Function)(void*);
}

#if 0
class ThreadBase;
#endif
class Thread;
class FunctionThread;
class Mutex;
class ConditionVariable;
class Barrier;
class SingleBarrier;
class Semaphore;
class SemaphoreB;
class RecursiveMutex;
class RWLock;
class Gate;

template <class M> class Lock;
template <class M> class TryLock;

template <> class Lock<Semaphore>;

template <class T, class M> class SafeVariable;
template <class T> class ThreadSpecificData;
template <> class ThreadSpecificData<void>;


//
// Mutex
//

//
// This Mutex class is designed to be used in a very stylized manner.
// While we can instantiate a Mutex, to lock or unlock it we must instantiate
// a Lock object with the Mutex.  This makes locked Mutexes exception-safe.
// If we want to trylock it we must use a TryLock object.
//

class Mutex
{
public:
    class Attr;
    Mutex();
    explicit Mutex(const Attr& attr);
    virtual ~Mutex();
    void lock();
    void unlock();
    bool trylock();
    pthread_mutex_t& theMutex();
private:
    Mutex(const Mutex&);
    void operator=(const Mutex&);
    void* operator new(size_t size);
    pthread_mutex_t m_mutex;
};


//
// Lock
//

//
// In order to ensure that Mutexes get released when an exception happens,
// we introduce yet another class that must be used to lock and unlock a
// Mutex.
//
// We have the following stylized method of defining a Mutex and then
// locking and unlocking it:
//
//   int   counter;       // A shared counter.
//   Mutex counterMT;     // A mutex on the counter.
//
//   {
//       //
//       // Instantiate Lock object in its own scope.
//       // This locks the Mutex `counterMT'.
//       //
//       Lock lock(counterMT);
//       //
//       // Update the counter.
//       //
//       counter += 1;
//       //
//       // On exit from the scope the Lock destructor is called which
//       // will unlock the Mutex `counterMT'.
//   }
//

template <class M = Mutex>
class Lock
{
public:
    explicit Lock(M& mutex);
    ~Lock();
private:
    Lock(const Lock&);
    void operator=(const Lock&);
    void* operator new(size_t size);
    M& m_mutex;
};

template <class M>
Lock<M>::Lock (M& mutex)
    : m_mutex(mutex)
{
    m_mutex.lock();
}

template <class M>
Lock<M>::~Lock()
{
    m_mutex.unlock();
}


//
// TryLock
//

//
// We have the following stylized method of defining a Mutex and then
// attempting to lock and unlock it:
//
//   int   counter;       // A shared counter.
//   Mutex counterMT;     // A mutex on the counter.
//
//   {
//       //
//       // Instantiate TryLock object in its own scope.
//       // This attempts to lock the Mutex `counterMT'.
//       //
//       TryLock trylock(counterMT);
//
//       if ( trylock.locked() )
//       {
//           //
//           // Update the counter.
//           //
//           counter += 1;
//       }
//       //
//       // On exit from the scope the TryLock destructor is called which
//       // will unlock the Mutex `counterMT' only if it was successfully
//       // locked.
//       //
//   }
//

template <class M = Mutex>
class TryLock
{
public:
    explicit TryLock(M& mutex);
    ~TryLock();
    bool locked() const;
private:
    TryLock(const TryLock&);
    void operator=(const TryLock&);
    void* operator new(size_t size);
    M& m_mutex;
    bool m_locked;
};


template <class M>
TryLock<M>::TryLock(M& mutex)
    : m_mutex(mutex), m_locked(mutex.trylock())
{
}

template <class M>
TryLock<M>::~TryLock()
{
    if ( m_locked )
    {
	m_mutex.unlock();
    }
}

template <class M>
bool
TryLock<M>::locked() const
{
    return m_locked;
}


//
// ConditionVariable
//

class ConditionVariable
    : public Mutex
{
public:
    class Attr;
    ConditionVariable();
    explicit ConditionVariable(const Attr& attr);
    ~ConditionVariable();
    void signal();		// Signal at least one waiting thread.
    void broadcast();		// Signal all waiting threads.
    void wait();			// Wait on the condition.
    void wait(Mutex& m);			// Wait on the condition.
    bool timedwait(const BoxLib::Time& abstime); // Wait until a specified time.
    bool timedwait(Mutex& m, const BoxLib::Time& abstime); // Wait until a specified time.
    pthread_cond_t& theCV();
private:
    ConditionVariable(const ConditionVariable&);
    void operator=(const ConditionVariable&);
    pthread_cond_t m_cv;
};


//
// Barrier
//

class Barrier
    : protected ConditionVariable
{
public:
    explicit Barrier(int cnt = 0);
    void init(int cnt);
    void wait();
private:
    Barrier(const Barrier&);
    void operator=(const Barrier&);
    int count;			// Number of threads to wait for
    int n_sleepers;		// Number of threads to waiting
    bool releasing;		// Still waking up sleepers
};


//
// Semaphore
//

class Semaphore
    : protected ConditionVariable
{
public:
    explicit Semaphore(int val_ = 1);
    void wait();
    bool trywait();
    void post();
private:
    Semaphore(const Semaphore&);
    void operator=(const Semaphore&);
    int value;
};

template <>
class Lock<Semaphore>
{
public:
    explicit Lock(Semaphore& sem_);
    ~Lock();
private:
    Semaphore& sem;
};


//
// SemaphoreB
//

class SemaphoreB
    : protected ConditionVariable
{
public:
    explicit SemaphoreB(int val_ = 1);
    int down();
    int up();
    int value();
    int decrement();
private:
    SemaphoreB(const SemaphoreB&);
    void operator=(const SemaphoreB&);
    int val;
};


//
//Single Barrier
//

class SingleBarrier
    : protected ConditionVariable
{
public:
    explicit SingleBarrier(int);
    void wait();
    void post();
private:
    SingleBarrier(const SingleBarrier&);
    void operator=(const SingleBarrier&);
    int count;			// Number of threads to post
    int n_posters;		// Number of threads who posted
    int n_waiters;		// Number of threads waiting
    bool releasing;		// Still waking up sleepers
};


//
//Recursive Mutex
//

class RecursiveMutex
    : protected ConditionVariable
{
public:
    RecursiveMutex();
    ~RecursiveMutex();
    void lock();
    bool trylock();
    void unlock();
private:
    RecursiveMutex(const RecursiveMutex&);
    void operator=(const RecursiveMutex&);
    int count;			// times locked recursively
    bool owned;
    pthread_t m_tid;
};


//
// Gate
//

class Gate
    : protected ConditionVariable
{
public:
    Gate();
    void open();
    void close();
    void release();
    void wait();
private:
    Gate(const Gate&);
    void operator=(const Gate&);
    bool closed;
};


//
// RWLock
//

class RWLock
{
public:
    RWLock();
    ~RWLock();
    void rdlock();
    bool tryrdlock();
    void wrlock();
    bool trywrlock();
    void unlock();
private:
    static void _wrlock_cleanup(void*);
    static void _rdlock_cleanup(void*);
    RWLock(const RWLock&);
    void operator=(const RWLock&);
    ConditionVariable readers;
    ConditionVariable writers;
    int state;			// -1:writer,0:free,>0:readers
    int waiters;			// number of waiting writers
};


//
//SafeVariable
//

template <class T, class M = Mutex>
class SafeVariable
{
public:
    void set(const T& n);
    const T get() const;

    SafeVariable& operator=(const T& n);
    operator const T() const;
private:
    T val;
    mutable M m;
};

template <class T, class M>
void
SafeVariable<T,M>::set(const T& n)
{
    m.lock();
    val = n;
    m.unlock();
}

template < class T, class M>
const T
SafeVariable<T,M>::get() const
{
    m.lock();
    const T tval = val;
    m.unlock();
    return tval;
}

template <class T, class M>
SafeVariable<T,M>&
SafeVariable<T,M>::operator=(const T& n)
{
    set(n);
    return *this;
}

template <class T, class M>
SafeVariable<T,M>::operator const T() const
{
    return get();
}


//
//ThreadSpecific Data
//

template <>
class ThreadSpecificData<void>
{
public:
    explicit ThreadSpecificData(void (*dst)(void*) = 0);
    virtual ~ThreadSpecificData() = 0;
    void* set(const void* t);
    void* get() const;
private:
    pthread_key_t key;
#ifndef BL_PTHREADS
    void* v;
    void (*tsd)(void*);
#endif
};

template <class T>
class ThreadSpecificData
    : private ThreadSpecificData<void>
{
public:
    explicit ThreadSpecificData(const T* p = 0);
    explicit ThreadSpecificData(const T* p, void (*dst)(void*) = 0);
    virtual ~ThreadSpecificData();
    T* set(const T* t);
    T* get() const;
    T* operator->();
    T& operator*();
    T* release();
    void reset(const T* p = 0);
private:
    ThreadSpecificData(const ThreadSpecificData&);
    void operator=(const ThreadSpecificData&);
    static void tsd_destroy(void*);
};

template <class T>
ThreadSpecificData<T>::ThreadSpecificData(const T* p)
    : ThreadSpecificData<void>(tsd_destroy)
{
    if ( p ) set(p);
}

template <class T>
ThreadSpecificData<T>::ThreadSpecificData(const T* p, void (*dst)(void*))
    : ThreadSpecificData<void>(dst)
{
    if ( p ) set(p);
}

template <class T>
ThreadSpecificData<T>::~ThreadSpecificData()
{
    reset();
}

template <class T>
T*
ThreadSpecificData<T>::set(const T* v)
{
    return static_cast<T*>(ThreadSpecificData<void>::set(static_cast<const void*>(v)));
}

template <class T>
T*
ThreadSpecificData<T>::get() const
{
    return static_cast<T*>(ThreadSpecificData<void>::get());
}

template <class T>
T*
ThreadSpecificData<T>::operator->()
{
    return get();
}

template <class T>
T&
ThreadSpecificData<T>::operator*()
{
    return *get();
}

template <class T>
T*
ThreadSpecificData<T>::release()
{
    return set(0);
}

template <class T>
void
ThreadSpecificData<T>::reset(const T* p)
{
    delete set(p);
}

template <class T>
void
ThreadSpecificData<T>::tsd_destroy(void* v)
{
    delete static_cast<T*>(v);
}

//
// Thread
//

class Thread
{
public:
    class thread_error;
    class Attr;
    //
    // Creates a Pthread.  Thread must be created, in the non-running state
    // because of the way constructors are implemented: virtual functions
    // can't be called in the constructor.
    //
    Thread();
    //
    // The destructor.
    //
    // This detach()s the Thread if it has neither been join()d or detach()d.
    // It doesn't cancel the thread.
    //
    virtual ~Thread();
    //
    // Start the thread running.
    //
    // Creates the POSIX thread and effectively executes this->work()
    // as its thread start routine.  It is an error to call this routine
    // more than once.  Note that this must be the first function called
    // on a Thread.
    //
    enum DetachState { Joinable, Detached };
    void run(DetachState instate = Joinable);
    void run(const Attr& attr);
    //
    // Detach ourselves.  A thread must either be detach()d or it
    // must be join()d by another thread to prevent memory leaks
    // in the POSIX thread layer.  This can only be called after
    // run() has been called.  It is an error to detach() a Thread
    // more than once or to detach() a Thread that has already
    // been join()d.
    //
    void detach();
    //
    // Detach specific thread.  A thread must either be detach()d or it
    // must be join()d by another thread to prevent memory leaks
    // in the POSIX thread layer.  This can only be called after
    // run() has been called on `thrd'.  It is an error to detach() a
    // Thread more than once or to detach() a Thread that has already
    // been join()d.
    //
    //
    // Go to sleep for the specified Time
    //
    static void sleep(const BoxLib::Time& tspec);
    //
    //
    //
    // Wait for a specific thread to terminate returning the exit status.
    // A thread must either be detach()d or it must be join()d by another
    // thread to prevent memory leaks in the POSIX thread layer.
    // It is an error to join() a Thread more than once or to join() a
    // Thread that has already been detach()d.
    //
    void* join() const;
    //
    // Attempt to cancel the specified thread; i.e. request that the
    // specified thread terminate gracefully.
    //
    void cancel();

    int getID() const;

    static int max_threads();
    //
    // Calling thread exits with passed return status.
    //
    static void exit(void* status = 0);
    //
    //
    // The calling thread will yield the CPU.
    //
    static void yield();
    //
    // Set the calling threads cancel state returning the previously set
    // state.  Valid values for `state' are Enable and Disable.
    //
    enum CancelState { Enable, Disable };
    static CancelState setCancelState(CancelState state);
    //
    // Set the calling threads cancel type returning the previously set
    // type.  Valid values for `type' are Deferred and Asynchronous.
    //
    enum CancelType { Deferred, Asynchronous };
    static CancelType setCancelType(CancelType type);
    //
    // Creates a deferred cancellation point in the calling thread.
    //
    static void testcancel();
    //
    // The thread-start routine that must be over-ridden by derived classes.
    // The return value is the exit value of the POSIX thread.
    //
    pthread_t& theThread();
    //
    // Self returns 0 if this is the 'master' thread.
    //
    static Thread* self();
protected:
    //
    // The following protected members should only be called by the executing thread.
    //
    virtual void* work() = 0;
private:
    static void* _doit(void* arg);
    //
    // Used by run() -- actual thread-start routine passed to pthread_create().
    //
    //    static void* DoIt (void* arg);
    //
    // Disallow copy constructor and assignment operator.
    //
    Thread(const Thread&);
    void operator=(const Thread&);
    static int get_next_id();
private:
    pthread_t m_tid;
    enum Status { Created, Running, Stopped };
    Status m_status;
    // a unique number on this processor, equal n'th thread that has started
    static int m_next_id;
    int m_id;
    mutable bool m_jod;
    static Mutex next_id;
    Mutex m_mutex;
    static ThreadSpecificData<Thread> m_self;
};

#if 0
class ThreadBase
{
protected:
    pthread_t m_tid;
    enum Status { Created, Running, Stopped };
    Status m_status;
    // a unique number on this processor, equal n'th thread that has started
    static int m_next_id;
    int m_id;
    mutable bool m_jod;
    static Mutex next_id;
    Mutex m_mutex;
    static ThreadSpecificData<ThreadBase> m_self;
private:
    static int get_next_id();
    ThreadBase(const ThreadBase&);
    void operator=(const ThreadBase&);
};
#endif

class Thread::thread_error
    : public std::runtime_error
{
public:
    thread_error(const char* file, int line, const char* call, int status = 0);
};



class FunctionThread
{
public:
    FunctionThread(Thread_Function func_, void* arg_ = 0, Thread::DetachState = Thread::Joinable);
    ~FunctionThread();
    void* join() const;
    void detach() const;
protected:
    mutable bool m_jod;
    pthread_t m_tid;
};

namespace BoxLib
{
    void Thread_Error(const char* file, int lineno, const char* message, int status);
}

#ifdef BL_PTHREADS

#define THREAD_REQUIRE(x)						\
do									\
{									\
  if ( int status = (x) )						\
    {									\
      BoxLib::Thread_Error(__FILE__, __LINE__, #x, status); 		\
    }									\
}									\
while ( false )


//
// Thread Attribute
//

class Thread::Attr
{
    friend class Thread;
public:
    explicit Attr(DetachState instate = Joinable);
    ~Attr();

    bool isJoinable() const;
#ifdef _POSIX_THREAD_PRIORITY_SCHEDULING
    enum SchedPolicy{ Regular, RoundRobin, FiFo };
    SchedPolicy set_schedpolicy(SchedPolicy policy);
    SchedPolicy get_schedpolicy() const;

    enum InheritSched { Inherit, Explicit };
    InheritSched set_inheritsched(InheritSched inherit);
    InheritSched get_inheritsched() const;

    enum SchedScope { System, Process };
    SchedScope set_scope(SchedScope scope);
    SchedScope get_scope() const;

    sched_param set_schedparam(const sched_param& param);
    sched_param get_schedparam() const;
#endif
#ifdef _POSIX_THREAD_ATTR_STACKADDR
    void* get_stackaddr() const;
    void* set_stackaddr(void*& addr);
#endif
#ifdef _POSIX_THREAD_ATTR_STACKSIZE
    size_t get_stacksize() const;
    size_t set_stacksize(size_t& size);
#endif
    pthread_attr_t theAttribute() const;
private:
    bool m_is_joinable;
    DetachState set_detachstate(DetachState detachstate);
    DetachState get_detachstate() const;

    Attr(const Attr&);
    void operator=(const Attr&);
    pthread_attr_t m_attr;
};

inline
pthread_attr_t
Thread::Attr::theAttribute() const
{
    return m_attr;
}
//
// Mutex Attributes
//

class Mutex::Attr
{
    friend class Mutex;
public:
    Attr();
    ~Attr();
    const pthread_mutexattr_t& attribute() const;
    pthread_mutexattr_t& attribute();
#ifdef _POSIX_THREAD_PROCESS_SHARED
    enum ProcessShared { Shared, Private };
    ProcessShared set_pshared(ProcessShared v);
    ProcessShared get_pshared() const;
#endif
private:
    Attr(const Attr&);
    void operator=(const Attr&);
    pthread_mutexattr_t m_attr;
};

//
// ConditionVariable Attributes
//

class ConditionVariable::Attr
{
    friend class ConditionVariable;
public:
    //
    // The default constructor.
    //
    Attr();
    //
    // The destructor.
    //
    ~Attr();
    //
    // Returns a reference to the pthread_condattr_t that can be
    // used in various portable and non-portable POSIX 1003.1c calls
    // that set condition variable attribute values.
    //
    const pthread_condattr_t& attribute() const;
    pthread_condattr_t& attribute();

#ifdef _POSIX_THREAD_PROCESS_SHARED
    enum ProcessShared { Shared, Private };
    ProcessShared set_pshared(ProcessShared v);
    ProcessShared get_pshared() const;
#endif

private:
    Attr(const Attr&);
    void operator=(const Attr&);
    pthread_condattr_t m_attr;
};

#else

#include <cstdlib>
//
// Thread Attribute
//

class Thread::Attr
{
    friend class Thread;
public:
    explicit Attr(DetachState instate = Joinable);
    bool isJoinable() const;
private:
    bool m_is_joinable;
    DetachState set_detachstate(DetachState detachstate);
    DetachState get_detachstate() const;
    Attr(const Attr&);
    void operator=(const Attr&);
};

//
// Mutex Attributes
//

class Mutex::Attr
{
    friend class Mutex;
public:
    Attr() {}

private:
    Attr(const Attr&);
    void operator=(const Attr&);
};

//
// ConditionVariable Attributes
//

class ConditionVariable::Attr
{
    friend class ConditionVariable;
public:
    Attr(){}
private:
    Attr(const Attr&);
    void operator=(const Attr&);
};


// Thread
inline Thread::~Thread() {}
inline int Thread::max_threads() { return 1;}
inline void Thread::yield() {}
inline void Thread::testcancel() {}
inline void* Thread::join() const { return 0;}

// Mutex
inline Mutex::Mutex() {}
inline Mutex::~Mutex() {}
inline void Mutex::lock() {}
inline void Mutex::unlock() {}

// ConditionVariable
inline ConditionVariable::ConditionVariable() {}
inline ConditionVariable::~ConditionVariable() {}
inline void ConditionVariable::broadcast() {}
inline void ConditionVariable::signal() {}
inline void ConditionVariable::wait() {}
inline void ConditionVariable::wait(Mutex&) {}

#endif

#endif
