#ifndef AMREX_INTERP_1D_C_H_
#define AMREX_INTERP_1D_C_H_
#include <AMReX_Config.H>

#include <AMReX_FArrayBox.H>
#include <AMReX_BCRec.H>
#include <AMReX_Vector.H>
#include <AMReX_Array.H>
#include <cmath>

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
pcinterp_interp (Box const& bx,
                 Array4<Real> const& fine, const int fcomp, const int ncomp,
                 Array4<Real const> const& crse, const int ccomp,
                 IntVect const& ratio) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    for (int n = 0; n < ncomp; ++n) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            const int ic = amrex::coarsen(i,ratio[0]);
            fine(i,0,0,n+fcomp) = crse(ic,0,0,n+ccomp);
        }
    }
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
nodebilin_slopes (Box const& bx, Array4<T> const& slope, Array4<T const> const& u,
                  const int icomp, const int ncomp, IntVect const& ratio) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    const Real rx = Real(1.)/ratio[0];

    for (int n = 0; n < ncomp; ++n) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            T dx0 = u(i+1,0,0,n+icomp) - u(i,0,0,n+icomp);
            slope(i,0,0,n) = rx*dx0;
        }
    }
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
nodebilin_interp (Box const& bx, Array4<T> const& fine, const int fcomp, const int ncomp,
                  Array4<T const> const& slope, Array4<T const> const& crse,
                  const int ccomp, IntVect const& ratio) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);
    const auto chi = amrex::ubound(slope);

    for (int n = 0; n < ncomp; ++n) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            const int ic = amrex::min(amrex::coarsen(i,ratio[0]),chi.x);
            const Real fx = i - ic*ratio[0];
            fine(i,0,0,n+fcomp) = crse(ic,0,0,n+ccomp) + fx*slope(ic,0,0,0);
        }
    }
}

// Or, compile time error this?
// Remove these functions in a way that gives a good error message?

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
facediv_face_interp (int /*ci*/, int /*cj*/, int /*ck*/,
                     int /*nc*/, int /*nf*/, int /*idir*/,
                     Array4<T const> const& /*crse*/,
                     Array4<T> const& /*fine*/,
                     Array4<const int> const& /*mask*/,
                     IntVect const& /*ratio*/) noexcept
{
    amrex::Abort("No 1D version of FaceDiv exists.\n");
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
facediv_int (int /*ci*/, int /*cj*/, int /*ck*/, int /*nf*/,
             GpuArray<Array4<T>, AMREX_SPACEDIM> const& /*fine*/,
             IntVect const& /*ratio*/,
             GpuArray<Real, AMREX_SPACEDIM> const& /*cellSize*/) noexcept
{
    amrex::Abort("No 1D version of FaceDiv exists.\n");
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
face_linear_interp_x (int i, int /*j*/, int /*k*/, int n, Array4<T> const& fine,
                      Array4<T const> const& crse, IntVect const& ratio) noexcept
{
    const int ii = amrex::coarsen(i,ratio[0]);
    if (i-ii*ratio[0] == 0) {
        fine(i,0,0,n) = crse(ii,0,0,n);
    } else {
        Real const w = static_cast<Real>(i-ii*ratio[0]) * (Real(1.)/ratio[0]);
        fine(i,0,0,n) = (Real(1.)-w) * crse(ii,0,0,n) + w * crse(ii+1,0,0,n);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void ccquartic_interp (int nvar,
                       Box const& fine_bx,
                       Box const& fine_bx2,
                       Box const& /* crse_bx */,
                       Box const& crse_bx2,
                       IntVect const& ratio,
                       Array4<T const> const& crse,
                       Array4<T>       const& fine ) noexcept

{
    const Array1D<Real, -2, 2> cL = { -0.01171875,  0.0859375, 0.5, -0.0859375, 0.01171875 };

    IntVect fblo  = fine_bx.lbound();
    IntVect fbhi  = fine_bx.ubound();
    IntVect fb2lo = fine_bx2.lbound();
    IntVect fb2hi = fine_bx2.ubound();
    // IntVect cblo  = crse_bx.lbound();
    // IntVect cbhi  = crse_bx.ubound();
    IntVect cb2lo = crse_bx2.lbound();
    IntVect cb2hi = crse_bx2.ubound();

    Array1D<Real, fb2lo.x, fb2hi.x> ftmp;

    for (int n = 0; n < nvar; ++n) {
        for (int i = cb2lo.x, i <= cb2hi.x; ++i) {
            int ii = 2*i;
            ftmp(ii  ) = 2.0 * ( cL(-2)*crse(i-2,0,0,n)
                               + cL(-1)*crse(i-1,0,0,n)
                               + cL( 0)*crse(i  ,0,0,n)
                               + cL( 1)*crse(i+1,0,0,n)
                               + cL( 2)*crse(i+2,0,0,n) );
            ftmp(ii+1) = 2.0 * crse(i,0,0,n) - ftmp(ii);
        } // i
        for (int ii = fblo.x; ii <= fbhi.x; ++ii) {
            fine(ii,0,0,n) = ftmp(ii);
        } // ii
    } // n
}

} // namespace amrex

#endif
