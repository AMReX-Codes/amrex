#ifndef AMREX_INTERP_BNDRYDATA_3D_K_H_
#define AMREX_INTERP_BNDRYDATA_3D_K_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_IArrayBox.H>

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void interpbndrydata_o1 (int i, int j, int k, int n,
                         Array4<Real> const& bdry, int nb,
                         Array4<Real const> const& crse, int nc, Dim3 const& r) noexcept
{
    int ic = amrex::coarsen(i,r.x);
    int jc = amrex::coarsen(j,r.y);
    int kc = amrex::coarsen(k,r.z);
    bdry(i,j,k,n+nb) = crse(ic,jc,kc,n+nc);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void interpbndrydata_x_o3 (int i, int j, int k, int n,
                           Array4<Real> const& bdry, int nb,
                           Array4<Real const> const& crse, int nc, Dim3 const& r,
                           Array4<int const> const& mask, int not_covered) noexcept
{
    int ic = amrex::coarsen(i,r.x);
    int jc = amrex::coarsen(j,r.y);
    int kc = amrex::coarsen(k,r.z);

    int lo = (mask(i,j-r.y,k) == not_covered) ? jc-1 : jc;
    int hi = (mask(i,j+r.y,k) == not_covered) ? jc+1 : jc;
    Real fac = (hi == lo+1) ? 1.0_rt : 0.5_rt;
    Real dy = fac*(crse(ic,hi,kc,n+nc)-crse(ic,lo,kc,n+nc));
    Real dy2 = (hi==lo+2) ? 0.5_rt*(crse(ic,jc+1,kc,n+nc) - 2._rt*crse(ic,jc,kc,n+nc) + crse(ic,jc-1,kc,n+nc)) : 0._rt;

    lo = (mask(i,j,k-r.z) == not_covered) ? kc-1 : kc;
    hi = (mask(i,j,k+r.z) == not_covered) ? kc+1 : kc;
    fac = (hi == lo+1) ? 1.0_rt : 0.5_rt;
    Real dz = fac*(crse(ic,jc,hi,n+nc)-crse(ic,jc,lo,n+nc));
    Real dz2 = (hi==lo+2) ? 0.5_rt*(crse(ic,jc,kc+1,n+nc) - 2._rt*crse(ic,jc,kc,n+nc) + crse(ic,jc,kc-1,n+nc)) : 0._rt;

    Real dyz = (mask(i,j-r.y,k-r.z) == not_covered and mask(i,j+r.y,k-r.z) == not_covered and
                mask(i,j-r.y,k+r.z) == not_covered and mask(i,j+r.y,k+r.z) == not_covered)
        ? 0.25_rt*(crse(ic,jc+1,kc+1,n+nc)-crse(ic,jc-1,kc+1,n+nc)+crse(ic,jc-1,kc-1,n+nc)-crse(ic,jc+1,kc-1,n+nc))
        : 0.0;

    Real y = -0.5_rt + (j-jc*r.y+0.5_rt)/r.y;
    Real z = -0.5_rt + (k-kc*r.z+0.5_rt)/r.z;
    bdry(i,j,k,n+nb) = crse(ic,jc,kc,n+nc) + y*dy + (y*y)*dy2 + z*dz + (z*z)*dz2 + y*z*dyz;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void interpbndrydata_y_o3 (int i, int j, int k, int n,
                           Array4<Real> const& bdry, int nb,
                           Array4<Real const> const& crse, int nc, Dim3 const& r,
                           Array4<int const> const& mask, int not_covered) noexcept
{
    int ic = amrex::coarsen(i,r.x);
    int jc = amrex::coarsen(j,r.y);
    int kc = amrex::coarsen(k,r.z);

    int lo = (mask(i-r.x,j,k) == not_covered) ? ic-1 : ic;
    int hi = (mask(i+r.x,j,k) == not_covered) ? ic+1 : ic;
    Real fac = (hi == lo+1) ? 1.0_rt : 0.5_rt;
    Real dx = fac*(crse(hi,jc,kc,n+nc)-crse(lo,jc,kc,n+nc));
    Real dx2 = (hi==lo+2) ? 0.5_rt*(crse(ic+1,jc,kc,n+nc) - 2._rt*crse(ic,jc,kc,n+nc) + crse(ic-1,jc,kc,n+nc)) : 0._rt;

    lo = (mask(i,j,k-r.z) == not_covered) ? kc-1 : kc;
    hi = (mask(i,j,k+r.z) == not_covered) ? kc+1 : kc;
    fac = (hi == lo+1) ? 1.0_rt : 0.5_rt;
    Real dz = fac*(crse(ic,jc,hi,n+nc)-crse(ic,jc,lo,n+nc));
    Real dz2 = (hi==lo+2) ? 0.5_rt*(crse(ic,jc,kc+1,n+nc) - 2._rt*crse(ic,jc,kc,n+nc) + crse(ic,jc,kc-1,n+nc)) : 0._rt;

    Real dxz = (mask(i-r.x,j,k-r.z) == not_covered and mask(i+r.x,j,k-r.z) == not_covered and
                mask(i-r.x,j,k+r.z) == not_covered and mask(i+r.x,j,k+r.z) == not_covered)
        ? 0.25_rt*(crse(ic+1,jc,kc+1,n+nc)-crse(ic-1,jc,kc+1,n+nc)+crse(ic-1,jc,kc-1,n+nc)-crse(ic+1,jc,kc-1,n+nc))
        : 0.0;


    Real x = -0.5_rt + (i-ic*r.x+0.5_rt)/r.x;
    Real z = -0.5_rt + (k-kc*r.z+0.5_rt)/r.z;
    bdry(i,j,k,n+nb) = crse(ic,jc,kc,n+nc) + x*dx + (x*x)*dx2 + z*dz + (z*z)*dz2 + x*z*dxz;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void interpbndrydata_z_o3 (int i, int j, int k, int n,
                           Array4<Real> const& bdry, int nb,
                           Array4<Real const> const& crse, int nc, Dim3 const& r,
                           Array4<int const> const& mask, int not_covered) noexcept
{
    int ic = amrex::coarsen(i,r.x);
    int jc = amrex::coarsen(j,r.y);
    int kc = amrex::coarsen(k,r.z);

    int lo = (mask(i-r.x,j,k) == not_covered) ? ic-1 : ic;
    int hi = (mask(i+r.x,j,k) == not_covered) ? ic+1 : ic;
    Real fac = (hi == lo+1) ? 1.0_rt : 0.5_rt;
    Real dx = fac*(crse(hi,jc,kc,n+nc)-crse(lo,jc,kc,n+nc));
    Real dx2 = (hi==lo+2) ? 0.5_rt*(crse(ic+1,jc,kc,n+nc) - 2._rt*crse(ic,jc,kc,n+nc) + crse(ic-1,jc,kc,n+nc)) : 0._rt;

    lo = (mask(i,j-r.y,k) == not_covered) ? jc-1 : jc;
    hi = (mask(i,j+r.y,k) == not_covered) ? jc+1 : jc;
    fac = (hi == lo+1) ? 1.0_rt : 0.5_rt;
    Real dy = fac*(crse(ic,hi,kc,n+nc)-crse(ic,lo,kc,n+nc));
    Real dy2 = (hi==lo+2) ? 0.5_rt*(crse(ic,jc+1,kc,n+nc) - 2._rt*crse(ic,jc,kc,n+nc) + crse(ic,jc-1,kc,n+nc)) : 0._rt;

    Real dxy = (mask(i-r.x,j-r.y,k) == not_covered and mask(i+r.x,j-r.y,k) == not_covered and
                mask(i-r.x,j+r.y,k) == not_covered and mask(i+r.x,j+r.y,k) == not_covered)
        ? 0.25_rt*(crse(ic+1,jc+1,kc,n+nc)-crse(ic-1,jc+1,kc,n+nc)+crse(ic-1,jc-1,kc,n+nc)-crse(ic+1,jc-1,kc,n+nc))
        : 0.0;

    Real x = -0.5_rt + (i-ic*r.x+0.5_rt)/r.x;
    Real y = -0.5_rt + (j-jc*r.y+0.5_rt)/r.y;
    bdry(i,j,k,n+nb) = crse(ic,jc,kc,n+nc) + x*dx + (x*x)*dx2 + y*dy + (y*y)*dy2 + x*y*dxy;
}

}
#endif
