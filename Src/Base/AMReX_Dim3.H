#ifndef AMREX_DIM3_H_
#define AMREX_DIM3_H_
#include <AMReX_Config.H>

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>

#include <iosfwd>

namespace amrex {

struct Dim3 { int x; int y; int z; };
struct XDim3 { Real x; Real y; Real z; };

std::ostream& operator<< (std::ostream& os, const Dim3& d);

}

#endif
