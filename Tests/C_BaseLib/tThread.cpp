#include <Thread.H>
#include <cassert>
#include <iostream>
#include <cstdlib>

namespace testing
{
void prodcons_main();
void philosophers_main();
}

namespace
{
Mutex rand_mutex;
    Mutex print_mutex;
int random_l()
{
    rand_mutex.lock();
    int i = rand();
    rand_mutex.unlock();
    return i;
}

//extern Mutex Debug::print_mutex;
#define PRINTMSG(x)							\
do									\
  {									\
    Lock<Mutex> lock(print_mutex);	\
    x;									\
  }									\
while ( false )


const int N_DINERS = 5;	// n philosophers sharing n chopsticks
Mutex chopsticks[N_DINERS];

// At most n philosophers are allowed into the room, others would
// have to wait at the door. This restriction demonstrates the use
// of condition variables.

ConditionVariable room_cond;
int room_occupancy = 0;
}

class philosopher* phillies[N_DINERS];

class philosopher
    : public Thread
{
public:
    philosopher(int id);
protected:
    void* work();
private:
    const int id;
};

philosopher::philosopher(int id_)
    : id(id_)
{
    assert(id >=0 && id < N_DINERS);
    run(Detached);
}

void*
philosopher::work()
{
    int l = id;
    int r = l+1;
    if ( r == N_DINERS )
    {
	r = 0;
    }
    if ( l & 1 )
    {
	int t = l;
	l = r;
	r = t;
    }

    PRINTMSG( std::cout << "Philosopher #" << id << " has entered the room." << std::endl );

    int count = random_l() % 10 + 1;
    while ( count-- )
    {
	chopsticks[l].lock();
	chopsticks[r].lock();
	PRINTMSG( std::cout << "Philosopher #" << id << " is eating spaghetti now." << std::endl );
	Thread::sleep(BoxLib::Time(random_l()%2,random_l()%1000000000));
	chopsticks[l].unlock();
	chopsticks[r].unlock();
	PRINTMSG( std::cout << "Philosopher #" << id << " is pondering about life." << std::endl );
	Thread::sleep(BoxLib::Time(random_l()%2,random_l()%1000000000));
    }

    room_cond.lock();
    room_occupancy--;
    phillies[id] = 0;
    room_cond.unlock();
    room_cond.signal();
    PRINTMSG( std::cout << "Philosopher #" << id << " has left the room (" << room_occupancy << " left)." << std::endl );
    return 0;
}

void
testing::philosophers_main()
{
    if ( Thread::max_threads() == 1 )
    {
	BoxLib::Error("not going to work");
    }
    room_cond.lock();
    for ( int i=0; i<N_DINERS; i++ )
    {
	phillies[i] = new philosopher(i);
    }

    room_occupancy = N_DINERS;

    for (;;)
    {
	while ( room_occupancy == N_DINERS )
	{
	    PRINTMSG( std::cout << "main thread about to block " << room_occupancy << std::endl );
	    room_cond.wait();
	}

	// Hm.. someone has left the room.

	room_cond.unlock();

	// Sleep for a while and then create a new philosopher

	PRINTMSG(std::cout << "main thread sleep" << std::endl);
	Thread::sleep(BoxLib::Time(1,200000000));
	PRINTMSG(std::cout << "main thread wake up" << std::endl);

	room_cond.lock();
	int i;
	for ( i=0; i<N_DINERS; i++ )
	{
	    if ( phillies[i] == 0 )
	    {
		break;
	    }
	}
	if ( i == N_DINERS )
	{
	    PRINTMSG( std::cout << "Contrary to what I was told, no one has left the room!!!!\n" );
	    PRINTMSG( std::cout << "I give up!!!" << std::endl );
	    BoxLib::Error("philosophers");
	}
	phillies[i] = new philosopher(i);
	room_occupancy++;
    }
}

namespace
{
ConditionVariable full;
ConditionVariable empty;
int empty_flag = 1;
const char* message;
const char* msgs[] = { "wibble", "wobble", "jelly", "plate" };
}

class producer
    : public Thread
{
public:
    producer(const char* name_) : name(name_) {}
    virtual void* work();
private:
    const char* name;
};

class consumer
    : public Thread
{
public:
    consumer(const char* name_) : name(name_) {}
    virtual void* work();
private:
    const char* name;
};

void
testing::prodcons_main()
{
    if ( Thread::max_threads() == 1 )
    {
	BoxLib::Error("not going to work");
    }
    PRINTMSG( std::cout << "main: creating producer1\n" );
    producer p1("producer1");
    p1.run();
    PRINTMSG( std::cout << "main: creating producer2\n" );
    producer p2("producer2");
    p2.run();

    PRINTMSG( std::cout << "main: creating consumer1\n" );
    consumer c1("consumer1");
    c1.run();
    PRINTMSG( std::cout << "main: creating consumer2\n" );
    consumer c2("consumer2");
    c2.run();

    PRINTMSG( std::cout << "main: creating consumer3\n" );

    consumer c3("consumer3");
    c3.work();
}

void*
consumer::work()
{
    for (;;)
    {
	full.lock();
	BoxLib::Time t = BoxLib::Time::get_time(); // 1/2 second from now
	t += 0.5;
	while ( empty_flag )
	{
	    PRINTMSG( std::cout << name << ": waiting for message\n" );
	    if ( full.timedwait(t) )
	    {
		PRINTMSG( std::cout << name << ": timed out, trying again\n" );
		t = BoxLib::Time::get_time();
		t += 0.5;
	    }
	    else if ( empty_flag )
	    {
		PRINTMSG( std::cout << name << ": woken but message already comsumed\n" );
	    }
	}

	PRINTMSG( std::cout << name << ": got message: '" << message << "'\n" );
	empty_flag = 1;
	empty.signal();
	full.unlock();
	Thread::sleep(BoxLib::Time(random_l() % 2, 1000000 * (random_l() % 1000)));
    }
}

void*
producer::work()
{
    for (;;)
    {
	full.lock();
	while ( !empty_flag )
	{
	    PRINTMSG( std::cout << name << ": having to wait for consumer\n" );
	    empty.wait(full);
	}
	message = msgs[random_l() % 4];
	empty_flag = 0;
	full.signal();
	PRINTMSG( std::cout << name << ": put message: '" << message << "'\n" );
	full.unlock();
	Thread::sleep(BoxLib::Time(random_l() % 2, 1000000 * (random_l() % 500) + 500) );
    }
}

int
main(int argc, char** argv)
{
    BoxLib::Initialize(argc, argv);
    using namespace testing;
    prodcons_main();
    philosophers_main();
    BoxLib::Finalize();
}
