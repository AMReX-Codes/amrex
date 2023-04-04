#ifndef AMREX_INTERP_BNDRYDATA_1D_K_H_
#define AMREX_INTERP_BNDRYDATA_1D_K_H_
#include <AMReX_Config.H>

#include <AMReX_Array4.H>

namespace amrex {

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void interpbndrydata_o1 (int i, int /*j*/, int /*k*/, int n,
                         Array4<T> const& bdry, int nb,
                         Array4<T const> const& crse, int nc,
                         Dim3 const& r) noexcept
{
    int ic = amrex::coarsen(i,r.x);
    bdry(i,0,0,n+nb) = crse(ic,0,0,n+nc);
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void interpbndrydata_x_o3 (int i, int /*j*/, int /*k*/, int n,
                           Array4<T> const& bdry, int nb,
                           Array4<T const> const& crse, int nc, Dim3 const& r,
                           Array4<int const> const& /*mask*/, int /*not_covered*/) noexcept
{
    int ic = amrex::coarsen(i,r.x);
    bdry(i,0,0,n+nb) = crse(ic,0,0,n+nc);
}

}
#endif
