#ifndef AMREX_MULITIFAB_UTIL_C_H_
#define AMREX_MULITIFAB_UTIL_C_H_

#include <AMReX_Gpu.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_IArrayBox.H>

namespace amrex {

AMREX_GPU_HOST_DEVICE
void amrex_int_to_real (Box const& bx, FArrayBox& rfab, IArrayBox const& ifab,
                        int rcomp, int icomp, int ncomp);

}

#endif

