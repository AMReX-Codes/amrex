#include <cassert>
#include <list>
#include <cstdio>

#include <WorkQueue.H>

namespace
{
Mutex print_mutex;
}

#if 1
#define DPRINTF(arg)							\
do									\
  {									\
    Lock<Mutex> lock(print_mutex);			\
    std::printf("tid(%ld): ", pthread_self());				\
    std::printf arg;							\
  }									\
while (false)
#else
#define DPRINTF(arg)
#endif

namespace
{
Mutex rand_mutex;
int random_l()
{
    Lock<Mutex> lock(rand_mutex);
    int i = rand();
    return i;
}
}

namespace
{
const int ITERATIONS = 25;

struct power_t
    : public WorkQueue::task
{
    power_t(int value_, int power_) : value(value_), power(power_)
    {
	assert(value >= 0);
	assert(power >= 0);
    }
    virtual void run();
    int value;
    int power;
};

struct engine_t
{
    explicit engine_t(const pthread_t& thr_) : thread_id(thr_), calls(1) {}
    pthread_t thread_id;
    int calls;
};

Mutex engine_list_mutex;
std::list<engine_t*> engine_list;
WorkQueue workq;
}

void
Engine_destructor(void* value_ptr)
{
    engine_t* engine = static_cast<engine_t*>(value_ptr);
    Lock<Mutex> lock(engine_list_mutex);
    engine_list.push_back(engine);
}

namespace
{
ThreadSpecificData<engine_t> engine_key(0, Engine_destructor);

void
power_t::run()
{
    engine_t* engine = engine_key.get();
    if ( engine == 0 )
    {
	engine = new engine_t(pthread_self());
	engine_key.set(engine);
    }
    else
    {
	engine->calls++;
    }
    int result = 1;
    DPRINTF(("Engine: computing %d^%d\n", value, power));
    for ( int count = 1; count <= power; count++ )
    {
	result *= value;
    }
    DPRINTF(("Engine: result %d\n", result));
}
}

extern "C"
void*
WorkQueue_routine(void*)
{
    for ( int count = 0; count < ITERATIONS; count++ )
    {
	power_t* element = new power_t(random_l() % 20, random_l() % 7);
	DPRINTF(("Request: %d^%d\n", element->value, element->power));
	workq.add(element);
	sleep(random_l() % 5);
    }
    return 0;
}

int
main(int argc, char** argv)
{
    BoxLib::Initialize(argc, argv);

    workq.max_threads(4);
    FunctionThread ft(WorkQueue_routine);
    WorkQueue_routine(0);
    ft.join();
    workq.add(0);

    workq.wait();
    workq.drain();
    sleep(1);			// FIXME: wait for descructors to all be called

    int count = 0;
    int calls = 0;
    for ( std::list<engine_t*>::const_iterator it = engine_list.begin(); it != engine_list.end(); ++it )
    {
	engine_t* engine = *it;
	count++;
	calls += engine->calls;
	std::printf("engine %d(%ld): %d calls\n", count, engine->thread_id, engine->calls);
	delete engine;
    }
    std::printf("%d engine threads processed %d calls\n", count, calls);
    BoxLib::Finalize();
}
