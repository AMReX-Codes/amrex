#include <cstdlib>
#include <ctime>


#include <Utility.H>
#include <Timer.H>

namespace
{
const long billion = 1000000000L;
}

Time::Time()
{
    tv_sec = 0;
    tv_nsec = 0;
}

Time::Time(long s, long n)
{
    BL_ASSERT(s >= 0);
    BL_ASSERT(n >= 0);
    BL_ASSERT(n < billion);
    tv_sec = s;
    tv_nsec = n;
    normalize();
}

Time::Time(double d)
{
    tv_sec = long(d);
    tv_nsec = long((d-tv_sec)*billion);
    normalize();
}

double
Time::as_double() const
{
    return tv_sec + tv_nsec/double(billion);
}

long
Time::as_long() const
{
    return tv_sec + tv_nsec/billion;
}

Time&
Time::operator+=(const Time& r)
{
    tv_sec += r.tv_sec;
    tv_nsec += r.tv_nsec;
    normalize();
    return *this;
}

Time
Time::operator+(const Time& r) const
{
    Time result(*this);
    return result+=r;
}

void
Time::normalize()
{
    if ( tv_nsec > billion )
    {
	tv_nsec -= billion;
	tv_sec += 1;
    }
}

Time
Time::get_time()
{
    return Time(BoxLib::wsecond());
}
