#ifndef _BL_THREAD_H_
#define _BL_THREAD_H_

//
// $Id: Thread.H,v 1.19 2002-04-08 20:51:36 car Exp $
//

//#include <cstdio>

#include <Utility.H>

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

class FunctionThread;
class Mutex;
class ConditionVariable;

class Barrier;
class SingleBarrier;
class Semaphore;
class SemaphoreB;
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
    class Implementation;
    Mutex();
    virtual ~Mutex();
    void lock();
    void unlock();
    bool trylock();
private:
    Mutex(const Mutex&);
    void operator=(const Mutex&);
    void* operator new(size_t size);
protected:
    Implementation* m_impl;
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
{
public:
    class Implementation;
    ConditionVariable();
    ~ConditionVariable();
    void signal();		// Signal at least one waiting thread.
    void broadcast();		// Signal all waiting threads.
    void wait();			// Wait on the condition.
    void lock();
    void unlock();
    bool trylock();
private:
    ConditionVariable(const ConditionVariable&);
    void operator=(const ConditionVariable&);
    Implementation* m_impl;
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
    class Implementation;
    explicit ThreadSpecificData(void (*dst)(void*) = 0);
    virtual ~ThreadSpecificData() = 0;
    void* set(const void* t);
    void* get() const;
private:
    Implementation* m_impl;
};

template <class T>
class ThreadSpecificData
    : private ThreadSpecificData<void>
{
public:
    explicit ThreadSpecificData(const T* p = 0);
    ThreadSpecificData (const T* p, void (*TSD_DESTROY)(void*));
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
ThreadSpecificData<T>::ThreadSpecificData(const T* p, void (*THR_DESTROY)(void*))
    : ThreadSpecificData<void>(THR_DESTROY)
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
    //printf("ThreadSpecificData<T>::tsd_destroy(%p)\n",v);
    delete static_cast<T*>(v);
}

//
// Thread
//

namespace Thread
{
    void sleep (const BoxLib::Time& tspec);
    unsigned long max_threads ();
    void exit (void* status = 0);
    void yield ();
    int getID ();
    bool baseThread ();
    enum CancelState { Enable, Disable };
    CancelState setCancelState (CancelState state);
}



class FunctionThread
{
public:
    enum DetachState { Joinable, Detached };
    class Implementation;
    FunctionThread(Thread_Function func_,
                   void* arg_ = 0,
                   DetachState = Joinable,
                   int stacksize = 0);
    ~FunctionThread();
    void* join() const;
    void detach() const;
protected:
    Implementation* m_impl;
};

namespace BoxLib
{
    void Thread_Error(const char* file, int lineno, const char* message, int status);
}

#ifndef BL_THREADS

#include <cstdlib>


// Thread
namespace Thread
{
inline unsigned long max_threads() { return 1UL;}
inline int getID () { return 0; }
inline bool baseThread () { return true; }
}

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
inline void ConditionVariable::lock() {}
inline void ConditionVariable::unlock() {}
 
#endif

#endif
