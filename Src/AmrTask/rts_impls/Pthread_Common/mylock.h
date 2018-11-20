#ifndef MYLOCK
#define MYLOCK

#include <pthread.h>

class MyLock
{
    private: 
	pthread_mutex_t _lock;

    public:
	MyLock(){
            pthread_mutex_init(&_lock, NULL);
	}
	~MyLock(){
	    pthread_mutex_destroy(&_lock);
	}
	void lock()
	{
	    pthread_mutex_lock(&_lock);
	}
	void unlock()
	{
	    pthread_mutex_unlock(&_lock);
	}
};
#endif
