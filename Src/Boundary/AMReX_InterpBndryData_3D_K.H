#ifndef AMREX_INTERP_BNDRYDATA_3D_K_H_
#define AMREX_INTERP_BNDRYDATA_3D_K_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_IArrayBox.H>

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void interpbndrydata_x_o1 (int islo, int ic, int jc, int kc, int n,
                           Array4<Real> const& bdry, int nb,
                           Array4<Real const> const& crse, int nc, Dim3 const& r) noexcept
{
    Real c = crse(ic,jc,kc,n+nc);
    const int i = ic*r.x + islo*(r.x-1);
    for     (int koff = 0; koff < r.z; ++koff) {
        int      k    =        kc*r.z  + koff;
        for (int joff = 0; joff < r.y; ++joff) {
            int  j    =        jc*r.y  + joff;
            bdry(i,j,k,n+nb) = c;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void interpbndrydata_y_o1 (int islo, int ic, int jc, int kc, int n,
                           Array4<Real> const& bdry, int nb,
                           Array4<Real const> const& crse, int nc, Dim3 const& r) noexcept
{
    Real c = crse(ic,jc,kc,n+nc);
    const int j = jc*r.y + islo*(r.y-1);
    for     (int koff = 0; koff < r.z; ++koff) {
        int      k    =        kc*r.z  + koff;
        for (int ioff = 0; ioff < r.x; ++ioff) {
            int  i    =        ic*r.x  + ioff;
            bdry(i,j,k,n+nb) = c;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void interpbndrydata_z_o1 (int islo, int ic, int jc, int kc, int n,
                           Array4<Real> const& bdry, int nb,
                           Array4<Real const> const& crse, int nc, Dim3 const& r) noexcept
{
    Real c = crse(ic,jc,kc,n+nc);
    const int k = kc*r.z + islo*(r.z-1);
    for     (int joff = 0; joff < r.y; ++joff) {
        int      j    =        jc*r.y  + joff;
        for (int ioff = 0; ioff < r.x; ++ioff) {
            int  i    =        ic*r.x  + ioff;
            bdry(i,j,k,n+nb) = c;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void interpbndrydata_x_o3 (int islo, int ic, int jc, int kc, int n,
                           Array4<Real> const& bdry, int nb,
                           Array4<Real const> const& crse, int nc, Dim3 const& r,
                           Array4<int const> const& mask, int not_covered) noexcept
{
    const int i = ic*r.x + islo*(r.x-1);
    int j = jc*r.y;
    int k = kc*r.z;

    int lo = (mask(i,j-1  ,k) == not_covered) ? jc-1 : jc;
    int hi = (mask(i,j+r.y,k) == not_covered) ? jc+1 : jc;
    Real fac = (hi == lo+1) ? 1.0 : 0.5;
    Real dy = fac*(crse(ic,hi,kc,n+nc)-crse(ic,lo,kc,n+nc));
    Real dy2 = (hi==lo+2) ? 0.5*(crse(ic,jc+1,kc,n+nc) - 2.*crse(ic,jc,kc,n+nc) + crse(ic,jc-1,kc,n+nc)) : 0.;

    lo = (mask(i,j,k-1  ) == not_covered) ? kc-1 : kc;
    hi = (mask(i,j,k+r.z) == not_covered) ? kc+1 : kc;
    fac = (hi == lo+1) ? 1.0 : 0.5;
    Real dz = fac*(crse(ic,jc,hi,n+nc)-crse(ic,jc,lo,n+nc));
    Real dz2 = (hi==lo+2) ? 0.5*(crse(ic,jc,kc+1,n+nc) - 2.*crse(ic,jc,kc,n+nc) + crse(ic,jc,kc-1,n+nc)) : 0.;

    Real dyz = (mask(i,j-1,k-1) == not_covered and mask(i,j+r.y,k-1) == not_covered
                and mask(i,j-1,k+r.z) == not_covered and mask(i,j+r.y,k+r.z) == not_covered)
        ? 0.25*(crse(ic,jc+1,kc+1,n+nc)-crse(ic,jc-1,kc+1,n+nc)+crse(ic,jc-1,kc-1,n+nc)-crse(ic,jc+1,kc-1,n+nc))
        : 0.0;

    Real c = crse(ic,jc,kc,n+nc);
    for     (int koff = 0; koff < r.z; ++koff) {
        k = kc*r.z  + koff;
        Real z = (0.5+(koff-r.z/2)) / r.z;
        for (int joff = 0; joff < r.y; ++joff) {
            j = jc*r.y  + joff;
            Real y = (0.5+(joff-r.y/2)) / r.y;
            bdry(i,j,k,n+nb) = c + y*dy + (y*y)*dy2 + z*dz + (z*z)*dz2 + y*z*dyz;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void interpbndrydata_y_o3 (int islo, int ic, int jc, int kc, int n,
                           Array4<Real> const& bdry, int nb,
                           Array4<Real const> const& crse, int nc, Dim3 const& r,
                           Array4<int const> const& mask, int not_covered) noexcept
{
    int i = ic*r.x;
    const int j = jc*r.y + islo*(r.y-1);
    int k = kc*r.z;

    int lo = (mask(i-1  ,j,k) == not_covered) ? ic-1 : ic;
    int hi = (mask(i+r.x,j,k) == not_covered) ? ic+1 : ic;
    Real fac = (hi == lo+1) ? 1.0 : 0.5;
    Real dx = fac*(crse(hi,jc,kc,n+nc)-crse(lo,jc,kc,n+nc));
    Real dx2 = (hi==lo+2) ? 0.5*(crse(ic+1,jc,kc,n+nc) - 2.*crse(ic,jc,kc,n+nc) + crse(ic-1,jc,kc,n+nc)) : 0.;

    lo = (mask(i,j,k-1  ) == not_covered) ? kc-1 : kc;
    hi = (mask(i,j,k+r.z) == not_covered) ? kc+1 : kc;
    fac = (hi == lo+1) ? 1.0 : 0.5;
    Real dz = fac*(crse(ic,jc,hi,n+nc)-crse(ic,jc,lo,n+nc));
    Real dz2 = (hi==lo+2) ? 0.5*(crse(ic,jc,kc+1,n+nc) - 2.*crse(ic,jc,kc,n+nc) + crse(ic,jc,kc-1,n+nc)) : 0.;

    Real dxz = (mask(i-1,j,k-1) == not_covered and mask(i+r.x,j,k-1) == not_covered
                and mask(i-1,j,k+r.z) == not_covered and mask(i+r.x,j,k+r.z) == not_covered)
        ? 0.25*(crse(ic+1,jc,kc+1,n+nc)-crse(ic-1,jc,kc+1,n+nc)+crse(ic-1,jc,kc-1,n+nc)-crse(ic+1,jc,kc-1,n+nc))
        : 0.0;

    Real c = crse(ic,jc,kc,n+nc);
    for     (int koff = 0; koff < r.z; ++koff) {
        k = kc*r.z  + koff;
        Real z = (0.5+(koff-r.z/2)) / r.z;
        for (int ioff = 0; ioff < r.x; ++ioff) {
            i = ic*r.x  + ioff;
            Real x = (0.5+(ioff-r.x/2)) / r.x;
            bdry(i,j,k,n+nb) = c + x*dx + (x*x)*dx2 + z*dz + (z*z)*dz2 + x*z*dxz;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void interpbndrydata_z_o3 (int islo, int ic, int jc, int kc, int n,
                           Array4<Real> const& bdry, int nb,
                           Array4<Real const> const& crse, int nc, Dim3 const& r,
                           Array4<int const> const& mask, int not_covered) noexcept
{
    int i = ic*r.x;
    int j = jc*r.y;
    const int k = kc*r.z + islo*(r.z-1);

    int lo = (mask(i-1  ,j,k) == not_covered) ? ic-1 : ic;
    int hi = (mask(i+r.x,j,k) == not_covered) ? ic+1 : ic;
    Real fac = (hi == lo+1) ? 1.0 : 0.5;
    Real dx = fac*(crse(hi,jc,kc,n+nc)-crse(lo,jc,kc,n+nc));
    Real dx2 = (hi==lo+2) ? 0.5*(crse(ic+1,jc,kc,n+nc) - 2.*crse(ic,jc,kc,n+nc) + crse(ic-1,jc,kc,n+nc)) : 0.;

    lo = (mask(i,j-1  ,k) == not_covered) ? jc-1 : jc;
    hi = (mask(i,j+r.y,k) == not_covered) ? jc+1 : jc;
    fac = (hi == lo+1) ? 1.0 : 0.5;
    Real dy = fac*(crse(ic,hi,kc,n+nc)-crse(ic,lo,kc,n+nc));
    Real dy2 = (hi==lo+2) ? 0.5*(crse(ic,jc+1,kc,n+nc) - 2.*crse(ic,jc,kc,n+nc) + crse(ic,jc-1,kc,n+nc)) : 0.;

    Real dxy = (mask(i-1,j-1,k) == not_covered and mask(i+r.x,j-1,k) == not_covered
                and mask(i-1,j+r.y,k) == not_covered and mask(i+r.x,j+r.y,k) == not_covered)
        ? 0.25*(crse(ic+1,jc+1,kc,n+nc)-crse(ic-1,jc+1,kc,n+nc)+crse(ic-1,jc-1,kc,n+nc)-crse(ic+1,jc-1,kc,n+nc))
        : 0.0;

    Real c = crse(ic,jc,kc,n+nc);
    for     (int joff = 0; joff < r.y; ++joff) {
        j = jc*r.y  + joff;
        Real y = (0.5+(joff-r.y/2)) / r.y;
        for (int ioff = 0; ioff < r.x; ++ioff) {
            i = ic*r.x + ioff;
            Real x = (0.5+(ioff-r.x/2)) / r.x;
            bdry(i,j,k,n+nb) = c + x*dx + (x*x)*dx2 + y*dy + (y*y)*dy2 + x*y*dxy;
        }
    }
}

}
#endif
