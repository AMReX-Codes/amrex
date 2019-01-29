#ifndef AMREX_INTERP_BNDRYDATA_1D_K_H_
#define AMREX_INTERP_BNDRYDATA_1D_K_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_IArrayBox.H>

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void interpbndrydata_x_o1 (int islo, int ic, int /*jc*/, int /*kc*/, int n, Array4<Real>& bdry,
                           Array4<Real> const& crse, Dim3 const& r)
{
    int i = ic*r.x + islo;
    bdry(i,0,0,n) = crse(ic,0,0,n);
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void interpbndrydata_x_o3 (int islo, int ic, int /*jc*/, int /*kc*/, int n, Array4<Real>& bdry,
                           Array4<Real> const& crse, Dim3 const& r,
                           Array4<int> const& mask, int not_covered)
{
    int i = ic*r.x + islo;
    bdry(i,0,0,n) = crse(ic,0,0,n);
}

}
#endif
