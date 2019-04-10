#ifndef AMREX_INTERP_BNDRYDATA_1D_K_H_
#define AMREX_INTERP_BNDRYDATA_1D_K_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_IArrayBox.H>

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void interpbndrydata_x_o1 (int islo, int ic, int /*jc*/, int /*kc*/, int n,
                           Array4<Real> const& bdry, int nb,
                           Array4<Real const> const& crse, int nc,
                           Dim3 const& r) noexcept
{
    const int i = ic*r.x + islo*(r.x-1);
    bdry(i,0,0,n+nb) = crse(ic,0,0,n+nc);
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void interpbndrydata_x_o3 (int islo, int ic, int /*jc*/, int /*kc*/, int n,
                           Array4<Real> const& bdry, int nb,
                           Array4<Real const> const& crse, int nc, Dim3 const& r,
                           Array4<int const> const& mask, int not_covered) noexcept
{
    const int i = ic*r.x + islo*(r.x-1);
    bdry(i,0,0,n+nb) = crse(ic,0,0,n+nc);
}

}
#endif
