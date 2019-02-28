#ifndef MD_K_H_
#define MD_K_H_

#include "Constants.H"

struct CheckPair
{
    template <class P>
    AMREX_GPU_DEVICE AMREX_INLINE
    bool operator()(const P& p1, const P& p2) const
    {
        amrex::Real d0 = (p1.pos(0) - p2.pos(0));
        amrex::Real d1 = (p1.pos(1) - p2.pos(1));
        amrex::Real d2 = (p1.pos(2) - p2.pos(2));    
        amrex::Real dsquared = d0*d0 + d1*d1 + d2*d2;   
        return (dsquared <= Params::cutoff*Params::cutoff);
    }
};

#endif
