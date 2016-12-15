#ifndef _PERIODICITY_H_
#define _PERIODICITY_H_

#include <IntVect.H>
#include <Box.H>
#include <vector>

// This provides length of period for periodic domains.  0 means it is 
// not periodic in that direction.
// It is also assumed that the periodic domain starts with index 0.
class Periodicity
{
public:
    Periodicity () : period(D_DECL(0,0,0)) {}
    explicit Periodicity (const IntVect& v) : period(v) {}

    bool isAnyPeriodic () const
	{ return D_TERM(period[0]>0, || period[1]>0, || period[2]>0); }
    bool isAllPeriodic () const
	{ return D_TERM(period[0]>0, && period[1]>0, && period[2]>0); }
    bool isPeriodic (int dir) const
	{ return period[dir]>0; }

    bool operator==(const Periodicity& rhs) const
	{ return period == rhs.period; }

    // Cell-centered domain Box "infinitely" long in non-periodic directions.
    Box Domain () const; 
    
    std::vector<IntVect> shiftIntVect () const;

    static const Periodicity& NonPeriodic ();

private:
    IntVect period;
};

#endif
