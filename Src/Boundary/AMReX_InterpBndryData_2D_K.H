#ifndef AMREX_INTERP_BNDRYDATA_2D_K_H_
#define AMREX_INTERP_BNDRYDATA_2D_K_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_IArrayBox.H>

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void interpbndrydata_x_o1 (int islo, int ic, int jc, int /*kc*/, int n, Array4<Real>& bdry,
                           Array4<Real> const& crse, Dim3 const& r)
{
    Real c = crse(ic,jc,0,n);
    int i = ic*r.x + islo;
    for (int joff = 0; joff < r.y; ++joff) {
        int  j    =        jc*r.y  + joff;
        bdry(i,j,0,n) = c;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void interpbndrydata_y_o1 (int islo, int ic, int jc, int /*kc*/, int n, Array4<Real>& bdry,
                           Array4<Real> const& crse, Dim3 const& r)
{
    Real c = crse(ic,jc,0,n);
    int j = jc*r.y + islo;
    for (int ioff = 0; ioff < r.x; ++ioff) {
        int  i    =        ic*r.x  + ioff;
        bdry(i,j,0,n) = c;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void interpbndrydata_x_o3 (int islo, int ic, int jc, int /*kc*/, int n, Array4<Real>& bdry,
                           Array4<Real> const& crse, Dim3 const& r,
                           Array4<int> const& mask, int not_covered)
{
    const int i = ic*r.x + islo;
    int j = jc*r.y;

    int lo = (mask(i,j-1  ,0) == not_covered) ? jc-1 : jc;
    int hi = (mask(i,j+r.y,0) == not_covered) ? jc+1 : jc;
    Real fac = (hi == lo+1) ? 1.0 : 0.5;
    Real dy = fac*(crse(ic,hi,0,n)-crse(ic,lo,0,n));
    Real dy2 = (hi==lo+2) ? 0.5*(crse(ic,jc+1,0,n) - 2.*crse(ic,jc,0,n) + crse(ic,jc-1,0,n)) : 0.;

    Real c = crse(ic,jc,0,n);
    for (int joff = 0; joff < r.y; ++joff) {
        j = jc*r.y  + joff;
        Real y = (0.5+(joff-r.y/2)) / r.y;
        bdry(i,j,0,n) = c + y*dy + (y*y)*dy2;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void interpbndrydata_y_o3 (int islo, int ic, int jc, int /*kc*/, int n, Array4<Real>& bdry,
                           Array4<Real> const& crse, Dim3 const& r,
                           Array4<int> const& mask, int not_covered)
{
    int i = ic*r.x;
    const int j = jc*r.y + islo;

    int lo = (mask(i-1  ,j,0) == not_covered) ? ic-1 : ic;
    int hi = (mask(i+r.x,j,0) == not_covered) ? ic+1 : ic;
    Real fac = (hi == lo+1) ? 1.0 : 0.5;
    Real dx = fac*(crse(hi,jc,0,n)-crse(lo,jc,0,n));
    Real dx2 = (hi==lo+2) ? 0.5*(crse(ic+1,jc,0,n) - 2.*crse(ic,jc,0,n) + crse(ic-1,jc,0,n)) : 0.;

    Real c = crse(ic,jc,0,n);
    for (int ioff = 0; ioff < r.x; ++ioff) {
        i = ic*r.x  + ioff;
        Real x = (0.5*(ioff-r.x/2)) / r.x;
        bdry(i,j,0,n) = c + x*dx + (x*x)*dx2;
    }
}

}
#endif
