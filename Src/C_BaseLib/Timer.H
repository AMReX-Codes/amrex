#ifndef _BL3_TIMER_H_
#define _BL3_TIMER_H_

#include <ctime>

#include <iostream>
#include <algorithm>

#include <Utility.H>

#ifdef WIN32
struct timespec
{
    long int tv_sec;
    long int tv_nsec;
};
#endif

class Time
    : public timespec
{
public:
    class bad_time;
    Time();
    Time(const timespec& ts);
    Time(long s, long n = 0);
    Time(double d);
    double as_double() const;
    long as_long() const;
    Time& operator+=(const Time&);
    Time operator+(const Time&) const;
    static Time get_time();
private:
    void normalize();
};

// Implements a simple in place Timer:
template <double (*FCN)(double*)> class base_Timer;

typedef base_Timer<BoxLib::wsecond> WallTimer;
typedef base_Timer<BoxLib::second>  CPUTimer;

template <double (*FCN)(double*)>
class base_Timer
{
public:
    class bad_timer;
    base_Timer();
    ~base_Timer();
    void start();
    void stop();
    void reset();
    double time() const;
    double accum_time() const;
    int count() const;
    bool is_running() const;
    static double tick();
private:
    bool running;
    double val;
    double held;
    double accum_val;
    int cnt;
};

template <double (*FCN)(double*)> std::ostream& operator<<(std::ostream&, const base_Timer<FCN>&);

template <double (*FCN)(double*)>
class base_Timer<FCN>::bad_timer
{
};

template <double (*FCN)(double*)>
base_Timer<FCN>::base_Timer()
    : running(false), val(0.0), accum_val(0.0), cnt(0)
{
}

template <double (*FCN)(double*)>
base_Timer<FCN>::~base_Timer()
{
    BLassert( !running );
}

template <double (*FCN)(double*)>
bool
base_Timer<FCN>::is_running() const
{
    return running;
}

template <double (*FCN)(double*)>
inline
void
base_Timer<FCN>::start()
{
    BL_ASSERT( !running );
    held = FCN(0);
    running = true;
}

template <double (*FCN)(double*)>
inline
void
base_Timer<FCN>::stop()
{
    BL_ASSERT( running );
    val = (FCN(0) - held);
#ifndef NDEBUG
    if ( val < 0 )
    {
	std::cout << "Got a negative time " << val << "!" << std::endl;
    }
#endif
    accum_val += val;
    cnt += 1;
    running = false;
}

template <double (*FCN)(double*)>
void
base_Timer<FCN>::reset()
{
    BL_ASSERT( !running );
    accum_val = 0;
    cnt = 0;
}

template <double (*FCN)(double*)>
double
base_Timer<FCN>::accum_time() const
{
    BL_ASSERT( !running );
    return accum_val;
}

template <double (*FCN)(double*)>
double
base_Timer<FCN>::time() const
{
    BL_ASSERT( !running );
    return val;
}

template <double (*FCN)(double*)>
int
base_Timer<FCN>::count() const
{
    BL_ASSERT( !running );
    return cnt;
}

// I stole this code from someplace, but I don't know where.
template <double (*FCN)(double*)>
double
base_Timer<FCN>::tick()
{
    const int M = 100;
    double timesfound[M];

    //  Collect a sequence of M unique time values from the system.

    for ( int i = 0; i < M; ++i )
    {
	double t2;
	double t1 = FCN(0);
	while( ((t2 = FCN(0)) - t1) == 0 )
	{
	}
	timesfound[i] = t2;
    }

    double minDelta = timesfound[1] - timesfound[0];
    for ( int i = 2; i < M; i++ )
    {
	minDelta = std::min(minDelta, std::max(timesfound[i]-timesfound[i-1], 0.0));
    }
    return minDelta;
}

template <double (*FCN)(double*)>
std::ostream& operator<<(std::ostream& os, const base_Timer<FCN>& bt)
{
    return os << "["
	      << bt.accum_time() << "/"
	      << bt.count()
	      << "]";
}

#endif
