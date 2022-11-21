#ifndef AMREX_INTERP_BNDRYDATA_3D_K_H_
#define AMREX_INTERP_BNDRYDATA_3D_K_H_
#include <AMReX_Config.H>

#include <AMReX_Array4.H>

namespace amrex {

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void interpbndrydata_o1 (int i, int j, int k, int n,
                         Array4<T> const& bdry, int nb,
                         Array4<T const> const& crse, int nc, Dim3 const& r) noexcept
{
    int ic = amrex::coarsen(i,r.x);
    int jc = amrex::coarsen(j,r.y);
    int kc = amrex::coarsen(k,r.z);
    bdry(i,j,k,n+nb) = crse(ic,jc,kc,n+nc);
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void interpbndrydata_x_o3 (int i, int j, int k, int n,
                           Array4<T> const& bdry, int nb,
                           Array4<T const> const& crse, int nc, Dim3 const& r,
                           Array4<int const> const& mask, int not_covered) noexcept
{
    int ic = amrex::coarsen(i,r.x);
    int jc = amrex::coarsen(j,r.y);
    int kc = amrex::coarsen(k,r.z);

    int lo = (mask(i,j-r.y,k) == not_covered) ? jc-1 : jc;
    int hi = (mask(i,j+r.y,k) == not_covered) ? jc+1 : jc;
    T fac = (hi == lo+1) ? T(1.0) : T(0.5);
    T dy = fac*(crse(ic,hi,kc,n+nc)-crse(ic,lo,kc,n+nc));
    T dy2 = (hi==lo+2) ? T(0.5)*(crse(ic,jc+1,kc,n+nc) - T(2.)*crse(ic,jc,kc,n+nc) + crse(ic,jc-1,kc,n+nc)) : T(0.);

    lo = (mask(i,j,k-r.z) == not_covered) ? kc-1 : kc;
    hi = (mask(i,j,k+r.z) == not_covered) ? kc+1 : kc;
    fac = (hi == lo+1) ? T(1.0) : T(0.5);
    T dz = fac*(crse(ic,jc,hi,n+nc)-crse(ic,jc,lo,n+nc));
    T dz2 = (hi==lo+2) ? T(0.5)*(crse(ic,jc,kc+1,n+nc) - T(2.)*crse(ic,jc,kc,n+nc) + crse(ic,jc,kc-1,n+nc)) : T(0.);

    T dyz = (mask(i,j-r.y,k-r.z) == not_covered && mask(i,j+r.y,k-r.z) == not_covered &&
                mask(i,j-r.y,k+r.z) == not_covered && mask(i,j+r.y,k+r.z) == not_covered)
        ? T(0.25)*(crse(ic,jc+1,kc+1,n+nc)-crse(ic,jc-1,kc+1,n+nc)+crse(ic,jc-1,kc-1,n+nc)-crse(ic,jc+1,kc-1,n+nc))
        : T(0.0);

    T y = -T(0.5) + (j-jc*r.y+T(0.5))/r.y;
    T z = -T(0.5) + (k-kc*r.z+T(0.5))/r.z;
    bdry(i,j,k,n+nb) = crse(ic,jc,kc,n+nc) + y*dy + (y*y)*dy2 + z*dz + (z*z)*dz2 + y*z*dyz;
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void interpbndrydata_y_o3 (int i, int j, int k, int n,
                           Array4<T> const& bdry, int nb,
                           Array4<T const> const& crse, int nc, Dim3 const& r,
                           Array4<int const> const& mask, int not_covered) noexcept
{
    int ic = amrex::coarsen(i,r.x);
    int jc = amrex::coarsen(j,r.y);
    int kc = amrex::coarsen(k,r.z);

    int lo = (mask(i-r.x,j,k) == not_covered) ? ic-1 : ic;
    int hi = (mask(i+r.x,j,k) == not_covered) ? ic+1 : ic;
    T fac = (hi == lo+1) ? T(1.0) : T(0.5);
    T dx = fac*(crse(hi,jc,kc,n+nc)-crse(lo,jc,kc,n+nc));
    T dx2 = (hi==lo+2) ? T(0.5)*(crse(ic+1,jc,kc,n+nc) - T(2.)*crse(ic,jc,kc,n+nc) + crse(ic-1,jc,kc,n+nc)) : T(0.);

    lo = (mask(i,j,k-r.z) == not_covered) ? kc-1 : kc;
    hi = (mask(i,j,k+r.z) == not_covered) ? kc+1 : kc;
    fac = (hi == lo+1) ? T(1.0) : T(0.5);
    T dz = fac*(crse(ic,jc,hi,n+nc)-crse(ic,jc,lo,n+nc));
    T dz2 = (hi==lo+2) ? T(0.5)*(crse(ic,jc,kc+1,n+nc) - T(2.)*crse(ic,jc,kc,n+nc) + crse(ic,jc,kc-1,n+nc)) : T(0.);

    T dxz = (mask(i-r.x,j,k-r.z) == not_covered && mask(i+r.x,j,k-r.z) == not_covered &&
                mask(i-r.x,j,k+r.z) == not_covered && mask(i+r.x,j,k+r.z) == not_covered)
        ? T(0.25)*(crse(ic+1,jc,kc+1,n+nc)-crse(ic-1,jc,kc+1,n+nc)+crse(ic-1,jc,kc-1,n+nc)-crse(ic+1,jc,kc-1,n+nc))
        : T(0.0);


    T x = -T(0.5) + (i-ic*r.x+T(0.5))/r.x;
    T z = -T(0.5) + (k-kc*r.z+T(0.5))/r.z;
    bdry(i,j,k,n+nb) = crse(ic,jc,kc,n+nc) + x*dx + (x*x)*dx2 + z*dz + (z*z)*dz2 + x*z*dxz;
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void interpbndrydata_z_o3 (int i, int j, int k, int n,
                           Array4<T> const& bdry, int nb,
                           Array4<T const> const& crse, int nc, Dim3 const& r,
                           Array4<int const> const& mask, int not_covered) noexcept
{
    int ic = amrex::coarsen(i,r.x);
    int jc = amrex::coarsen(j,r.y);
    int kc = amrex::coarsen(k,r.z);

    int lo = (mask(i-r.x,j,k) == not_covered) ? ic-1 : ic;
    int hi = (mask(i+r.x,j,k) == not_covered) ? ic+1 : ic;
    T fac = (hi == lo+1) ? T(1.0) : T(0.5);
    T dx = fac*(crse(hi,jc,kc,n+nc)-crse(lo,jc,kc,n+nc));
    T dx2 = (hi==lo+2) ? T(0.5)*(crse(ic+1,jc,kc,n+nc) - T(2.)*crse(ic,jc,kc,n+nc) + crse(ic-1,jc,kc,n+nc)) : T(0.);

    lo = (mask(i,j-r.y,k) == not_covered) ? jc-1 : jc;
    hi = (mask(i,j+r.y,k) == not_covered) ? jc+1 : jc;
    fac = (hi == lo+1) ? T(1.0) : T(0.5);
    T dy = fac*(crse(ic,hi,kc,n+nc)-crse(ic,lo,kc,n+nc));
    T dy2 = (hi==lo+2) ? T(0.5)*(crse(ic,jc+1,kc,n+nc) - T(2.)*crse(ic,jc,kc,n+nc) + crse(ic,jc-1,kc,n+nc)) : T(0.);

    T dxy = (mask(i-r.x,j-r.y,k) == not_covered && mask(i+r.x,j-r.y,k) == not_covered &&
                mask(i-r.x,j+r.y,k) == not_covered && mask(i+r.x,j+r.y,k) == not_covered)
        ? T(0.25)*(crse(ic+1,jc+1,kc,n+nc)-crse(ic-1,jc+1,kc,n+nc)+crse(ic-1,jc-1,kc,n+nc)-crse(ic+1,jc-1,kc,n+nc))
        : T(0.0);

    T x = -T(0.5) + (i-ic*r.x+T(0.5))/r.x;
    T y = -T(0.5) + (j-jc*r.y+T(0.5))/r.y;
    bdry(i,j,k,n+nb) = crse(ic,jc,kc,n+nc) + x*dx + (x*x)*dx2 + y*dy + (y*y)*dy2 + x*y*dxy;
}

}
#endif
