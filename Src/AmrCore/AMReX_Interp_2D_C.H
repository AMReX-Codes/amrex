#ifndef AMREX_INTERP_2D_C_H_
#define AMREX_INTERP_2D_C_H_
#include <AMReX_Config.H>

#include <AMReX_FArrayBox.H>
#include <AMReX_BCRec.H>
#include <AMReX_Vector.H>
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
        for (int j = lo.y; j <= hi.y; ++j) {
            const int jc = amrex::coarsen(j,ratio[1]);
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                const int ic = amrex::coarsen(i,ratio[0]);
                fine(i,j,0,n+fcomp) = crse(ic,jc,0,n+ccomp);
            }
        }
    }
}

namespace {
    static constexpr int ix   = 0;
    static constexpr int iy   = 1;
    static constexpr int ixy  = 2;
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
nodebilin_slopes (Box const& bx, Array4<T> const& slope, Array4<T const> const& u,
                  const int icomp, const int ncomp, IntVect const& ratio) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    const Real rx = Real(1.)/ratio[0];
    const Real ry = Real(1.)/ratio[1];

    for (int n = 0; n < ncomp; ++n) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                T dx0 = u(i+1,j,0,n+icomp) - u(i,j,0,n+icomp);
                T d0x = u(i,j+1,0,n+icomp) - u(i,j,0,n+icomp);
                T dx1 = u(i+1,j+1,0,n+icomp) - u(i,j+1,0,n+icomp);

                slope(i,j,0,n+ncomp*ix ) = rx*dx0;
                slope(i,j,0,n+ncomp*iy ) = ry*d0x;
                slope(i,j,0,n+ncomp*ixy) = rx*ry*(dx1 - dx0);
            }
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
        for (int j = lo.y; j <= hi.y; ++j) {
            const int jc = amrex::min(amrex::coarsen(j,ratio[1]),chi.y);
            const Real fy = j - jc*ratio[1];
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                const int ic = amrex::min(amrex::coarsen(i,ratio[0]),chi.x);
                const Real fx = i - ic*ratio[0];
                fine(i,j,0,n+fcomp) = crse(ic,jc,0,n+ccomp)
                    + fx*slope(ic,jc,0,n+ncomp*ix)
                    + fy*slope(ic,jc,0,n+ncomp*iy)
                    + fx*fy*slope(ic,jc,0,n+ncomp*ixy);
            }
        }
    }
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
facediv_face_interp (int ci, int cj, int /*ck*/,
                     int nc, int nf, int idir,
                     Array4<T const> const& crse,
                     Array4<T> const& fine,
                     Array4<int const> const& mask,
                     IntVect const& ratio) noexcept
{
    if (mask) {
        if (!mask(ci, cj, 0, nc))
            { return; }
    }

    const int fi = ci*ratio[0];
    const int fj = cj*ratio[1];

    switch (idir) {
        case 0:
        {
            const Real neg = crse(ci, cj-1, 0, nc);
            const Real cen = crse(ci, cj  , 0, nc);
            const Real pos = crse(ci, cj+1, 0, nc);

            fine(fi, fj  , 0, nf) = Real(0.125)*(8*cen + neg - pos);
            fine(fi, fj+1, 0, nf) = Real(0.125)*(8*cen + pos - neg);

            break;
        }
        case 1:
        {
            const Real neg = crse(ci-1, cj, 0, nc);
            const Real cen = crse(ci  , cj, 0, nc);
            const Real pos = crse(ci+1, cj, 0, nc);

            fine(fi  , fj, 0, nf) = Real(0.125)*(8*cen + neg - pos);
            fine(fi+1, fj, 0, nf) = Real(0.125)*(8*cen + pos - neg);

            break;
        }
    }
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
facediv_int (int ci, int cj, int /*ck*/, int nf,
             GpuArray<Array4<T>, AMREX_SPACEDIM> const& fine,
             IntVect const& ratio,
             GpuArray<Real, AMREX_SPACEDIM> const& cellSize) noexcept
{
    const int fi = ci*ratio[0];
    const int fj = cj*ratio[1];

    // References to fine exterior values.
    const Real umm = fine[0](fi,   fj,   0, nf);
    const Real ump = fine[0](fi,   fj+1, 0, nf);
    const Real upm = fine[0](fi+2, fj,   0, nf);
    const Real upp = fine[0](fi+2, fj+1, 0, nf);

    const Real vmm = fine[1](fi,   fj,   0, nf);
    const Real vmp = fine[1](fi+1, fj,   0, nf);
    const Real vpm = fine[1](fi,   fj+2, 0, nf);
    const Real vpp = fine[1](fi+1, fj+2, 0, nf);

    const Real dxdy = cellSize[0]/cellSize[1];
    const Real x_corr = Real(0.25)*dxdy * (vpp+vmm-vmp-vpm);
    const Real y_corr = Real(0.25)/dxdy * (upp+umm-ump-upm);

    // Calc fine faces on interior of coarse cells.
    fine[0](fi+1,fj  ,0,nf) = Real(0.5)*(umm+upm) + x_corr;
    fine[0](fi+1,fj+1,0,nf) = Real(0.5)*(ump+upp) + x_corr;
    fine[1](fi,  fj+1,0,nf) = Real(0.5)*(vmm+vpm) + y_corr;
    fine[1](fi+1,fj+1,0,nf) = Real(0.5)*(vmp+vpp) + y_corr;
}



template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
face_linear_interp_x (int i, int j, int /*k*/, int n, Array4<T> const& fine,
                      Array4<T const> const& crse, IntVect const& ratio) noexcept
{
    const int ii = amrex::coarsen(i,ratio[0]);
    const int jj = amrex::coarsen(j,ratio[1]);
    if (i-ii*ratio[0] == 0) {
        fine(i,j,0,n) = crse(ii,jj,0,n);
    } else {
        Real const w = static_cast<Real>(i-ii*ratio[0]) * (Real(1.)/ratio[0]);
        fine(i,j,0,n) = (Real(1.)-w) * crse(ii,jj,0,n) + w * crse(ii+1,jj,0,n);
    }
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
face_linear_interp_y (int i, int j, int /*k*/, int n, Array4<T> const& fine,
                      Array4<T const> const& crse, IntVect const& ratio) noexcept
{
    const int ii = amrex::coarsen(i,ratio[0]);
    const int jj = amrex::coarsen(j,ratio[1]);
    if (j-jj*ratio[1] == 0) {
        fine(i,j,0,n) = crse(ii,jj,0,n);
    } else {
        Real const w = static_cast<Real>(j-jj*ratio[1]) * (Real(1.)/ratio[1]);
        fine(i,j,0,n) = (Real(1.)-w) * crse(ii,jj,0,n) + w * crse(ii,jj+1,0,n);
    }
}

template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void ccprotect_check_redo (Box const& interp_bx,
                           bool& redo_me,
                           const int n,
                           Array4<T>       const& fine,
                           Array4<T const> const& fine_state) noexcept
{
    bool r = false;
    Dim3 lo = lbound(interp_bx), hi = ubound(interp_bx);
    for     (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
            r = r || ((fine_state(i,j,n) + fine(i,j,n)) < 0.0 );
        }
    }

    redo_me = r;
}

template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void ccprotect_copy_fine_orig(Box const& interp_bx,
                              const int n,
                              Array4<T> const& fine_orig,
                              Array4<T> const& fine) noexcept
{
    Dim3 lo = lbound(interp_bx), hi = ubound(interp_bx);
    for     (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
            fine_orig(i-lo.x,j-lo.y,0) = fine(i,j,n);
        }
    }
}

template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void ccprotect_calc_sums (Box const& interp_bx,
                          Real& crseTot,
                          Real& SumN,
                          Real& SumP,
                          const int n,
                          GpuArray<Vector<T>, 2> fvc,
                          Array4<T> const& fine,
                          Array4<T const> const& fine_state) noexcept
{
    Dim3 lo = lbound(interp_bx), hi = ubound(interp_bx);
    Real pos = 0.0;
    Real neg = 0.0;
    Real tot = 0.0;
    Real fvol;
    for     (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
            fvol = (fvc[0][i+1]-fvc[0][i]) * (fvc[1][j+1]-fvc[1][j]);
            tot += fvol * fine(i,j,n);
            if (fine_state(i,j,n) <= 0.0) {
                neg += fvol * fine_state(i,j,n);
            } else {
                pos += fvol * fine_state(i,j,n);
            }
        }
    }
    crseTot = tot;
    SumN = neg;
    SumP = pos;
}

template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void ccprotect_case1 (Box const& interp_bx,
                      const Real crseTot,
                      const Real SumN,
                      const Real SumP,
                      const int n,
                      const Real cvol,
                      Array4<T> const& fine,
                      Array4<T const> const& fine_state) noexcept
{
    Dim3 lo = lbound(interp_bx), hi = ubound(interp_bx);

    // Fill in negative values first.
    for     (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
            if (fine_state(i,j,n) < 0.0) {
                fine(i,j,n) = -fine_state(i,j,n);
            }
        }
    }

    // Then, add the remaining positive proportionally.
    if (SumP > 0.0) {
        Real alpha = (crseTot - Math::abs(SumN)) / SumP;
        for     (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                if (fine_state(i,j,n) > 0.0) {
                    fine(i,j,n) = alpha * fine_state(i,j,n);
                }
            }
        }
    } else { /* (SumP < 0) */
        Real posVal = (crseTot - Math::abs(SumN)) / cvol;
        for     (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                fine(i,j,n) = fine(i,j,n) + posVal;
            }
        }
    }
}

template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void ccprotect_case2 (Box const& interp_bx,
                      const Real crseTot,
                      const Real SumN,
                      const Real /* SumP */,
                      const int n,
                      Array4<T> const& fine,
                      Array4<T const> const& fine_state) noexcept
{
    Dim3 lo = lbound(interp_bx), hi = ubound(interp_bx);
    Real alpha = crseTot / Math::abs(SumN);
    for     (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
            if (fine_state(i,j,n) < 0.0) {
                // Add correction to negative cells proportionally.
                fine(i,j,n) = alpha * Math::abs(fine_state(i,j,n));
            } else {
                // Don't touch the positive states.
                fine(i,j,n) = 0.0;
            }
        }
    }
}

template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void ccprotect_case3 (Box const& interp_bx,
                      const Real crseTot,
                      const Real SumN,
                      const Real SumP,
                      const int n,
                      const Real cvol,
                      Array4<T> const& fine,
                      Array4<T const> const& fine_state) noexcept
{
    Dim3 lo = lbound(interp_bx), hi = ubound(interp_bx);
    Real negVal = (SumP + SumN + crseTot) / cvol;
    for     (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
            // We want to make them all as close to the same negative value.
            fine(i,j,n) = negVal - fine_state(i,j,n);
        }
    }
}

template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void ccprotect_case4 (Box const& interp_bx,
                      const Real crseTot,
                      const Real SumN,
                      const Real SumP,
                      const int n,
                      Array4<T> const& fine,
                      Array4<T const> const& fine_state) noexcept
{
    Dim3 lo = lbound(interp_bx), hi = ubound(interp_bx);
    Real alpha = (crseTot + SumN) / SumP;
    for     (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
            if (fine_state(i,j,n) < 0.0) {
                // Absorb the negative correction
                fine(i,j,n) = -fine_state(i,j,n);
            } else {
                // Redistribute the rest to make negative cells positive
                fine(i,j,n) = alpha * fine_state(i,j,n);
            }
        }
    }
}

}  // namespace amrex

#endif
