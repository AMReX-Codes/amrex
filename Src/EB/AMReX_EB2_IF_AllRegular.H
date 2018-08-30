#ifndef AMREX_EB2_IF_ALL_REGULAR_H_
#define AMREX_EB2_IF_ALL_REGULAR_H_

#include <AMReX_Array.H>

namespace amrex { namespace EB2 {

// For all implicit functions, >0: body; =0: boundary; <0: fluid

class AllRegularIF
{
public:
    constexpr Real operator() (const RealArray&) const { return -1.0; }
};

}}

#endif
