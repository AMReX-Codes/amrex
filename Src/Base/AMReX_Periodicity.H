#ifndef AMREX_PERIODICITY_H_
#define AMREX_PERIODICITY_H_

#include <AMReX_IntVect.H>
#include <AMReX_Box.H>
#include <vector>

namespace amrex {

/**
* \brief This provides length of period for periodic domains.  0 means it is 
* not periodic in that direction.
* It is also assumed that the periodic domain starts with index 0.
*/
class Periodicity
{
public:
    Periodicity () noexcept : period(AMREX_D_DECL(0,0,0)) {}
    explicit Periodicity (const IntVect& v) noexcept : period(v) {}

    bool isAnyPeriodic () const noexcept
	{ return AMREX_D_TERM(period[0]>0, || period[1]>0, || period[2]>0); }
    bool isAllPeriodic () const noexcept
	{ return AMREX_D_TERM(period[0]>0, && period[1]>0, && period[2]>0); }
    bool isPeriodic (int dir) const noexcept
	{ return period[dir]>0; }

    bool operator==(const Periodicity& rhs) const noexcept
	{ return period == rhs.period; }

    //! Cell-centered domain Box "infinitely" long in non-periodic directions.
    Box Domain () const noexcept;
    
    std::vector<IntVect> shiftIntVect () const;

    static const Periodicity& NonPeriodic () noexcept;

private:
    IntVect period;
};

}

#endif
