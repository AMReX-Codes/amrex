#ifndef AMREX_EB2_IF_ALL_REGULAR_H_
#define AMREX_EB2_IF_ALL_REGULAR_H_

#include <AMReX_Array.H>
#include <AMReX_EB2_IF_Base.H>

namespace amrex { namespace EB2 {

// For all implicit functions, >0: body; =0: boundary; <0: fluid

class AllRegularIF
    : public GPUable
{
public:
    constexpr Real operator() (const RealArray&) const noexcept { return -1.0; }

    AMREX_GPU_HOST_DEVICE
    constexpr Real operator() (AMREX_D_DECL(Real x, Real y, Real z)) const noexcept { return -1.0; }
};

}}

#endif
