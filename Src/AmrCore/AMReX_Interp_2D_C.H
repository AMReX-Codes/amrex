#ifndef AMREX_INTERP_2D_C_H_
#define AMREX_INTERP_2D_C_H_
#include <AMReX_Config.H>

#include <AMReX_Array.H>
#include <AMReX_Geometry.H>

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
    constexpr int ix   = 0;
    constexpr int iy   = 1;
    constexpr int ixy  = 2;
}

template<typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
nodebilin_slopes (Box const& bx, Array4<T> const& slope, Array4<T const> const& u,
                  const int icomp, const int ncomp, IntVect const& ratio) noexcept
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);

    const Real rx = Real(1.)/Real(ratio[0]);
    const Real ry = Real(1.)/Real(ratio[1]);

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

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void ccprotect_2d (int ic, int jc, int /*kc*/, int nvar,
                   Box const& fine_bx,
                   IntVect const& ratio,
                   GeometryData cs_geomdata,
                   GeometryData fn_geomdata,
                   Array4<T>       const& fine,
                   Array4<T const> const& fine_state) noexcept
{
    // Calculate bounds for interpolation
    Dim3 fnbxlo = lbound(fine_bx);
    Dim3 fnbxhi = ubound(fine_bx);
    int ilo = amrex::max(ratio[0]*ic,              fnbxlo.x);
    int ihi = amrex::min(ratio[0]*ic+(ratio[0]-1), fnbxhi.x);
    int jlo = amrex::max(ratio[1]*jc,              fnbxlo.y);
    int jhi = amrex::min(ratio[1]*jc+(ratio[1]-1), fnbxhi.y);

    /*
     * Check if interpolation needs to be redone for derived components (n > 0)
     */
    for (int n = 1; n < nvar-1; ++n) {

        bool redo_me = false;
        for     (int j = jlo; j <= jhi; ++j) {
            for (int i = ilo; i <= ihi; ++i) {
                if ((fine_state(i,j,0,n) + fine(i,j,0,n)) < 0.0) {
                    redo_me = true;
                }
            }
        }

        /*
         * If all the fine values are non-negative after the original
         * interpolated correction, then we do nothing here.
         *
         * If any of the fine values are negative after the original
         * interpolated correction, then we do our best.
         */
        if (redo_me) {

            /*
             * Calculate coarse cell volumes.
             */
            Real cvol;
            // Calculate coarse cell volume
            const Real* cs_dx = cs_geomdata.CellSize();
            if (cs_geomdata.coord == CoordSys::cartesian) {
                // Cartesian
                cvol = cs_dx[0] * cs_dx[1];
            } else {
                // RZ
                const Real* cs_problo = cs_geomdata.ProbLo();
                Real rp = cs_problo[0] + (ic+0.5_rt)*cs_dx[0];
                Real rm = cs_problo[0] + (ic-0.5_rt)*cs_dx[0];
                cvol = (rp*rp - rm*rm) * cs_dx[1];
            }

            /*
             * First, calculate the following quantities:
             *
             * crseTot = volume-weighted sum of all interpolated values
             *           of the correction, which is equivalent to
             *           the total volume-weighted coarse correction
             *
             * SumN = volume-weighted sum of all negative values of fine_state
             *
             * SumP = volume-weighted sum of all positive values of fine_state
             */
            Real crseTot = 0.0;
            Real SumN = 0.0;
            Real SumP = 0.0;
            for     (int j = jlo; j <= jhi; ++j) {
                for (int i = ilo; i <= ihi; ++i) {

                    Real fvol;
                    // Calculate fine cell volume
                    const Real* fn_dx = fn_geomdata.CellSize();
                    if (fn_geomdata.coord == CoordSys::cartesian) {
                        // Cartesian
                        fvol = fn_dx[0] * fn_dx[1];
                    } else {
                        // RZ
                        const Real* fn_problo = fn_geomdata.ProbLo();
                        Real rp = fn_problo[0] + (i+0.5_rt)*fn_dx[0];
                        Real rm = fn_problo[0] + (i-0.5_rt)*fn_dx[0];
                        fvol = (rp*rp - rm*rm) * fn_dx[1];
                    }

                    // Calculate sums
                    crseTot += fvol * fine(i,j,0,n);
                    if (fine_state(i,j,0,n) <= 0.0) {
                        SumN += fvol * fine_state(i,j,0,n);
                    } else {
                        SumP += fvol * fine_state(i,j,0,n);
                    }
                }
            }

            if ( (crseTot > 0) && (crseTot > std::abs(SumN)) ) {

                /*
                 * Special case 1:
                 *
                 * Coarse correction > 0, and fine_state has some cells
                 * with negative values which will be filled before
                 * adding to the other cells.
                 *
                 * Use the correction to bring negative cells to zero,
                 * then distribute the remaining positive proportionally.
                 */
                for     (int j = jlo; j <= jhi; ++j) {
                    for (int i = ilo; i <= ihi; ++i) {

                        // Fill in negative values first.
                        if (fine_state(i,j,0,n) < 0.0) {
                            fine(i,j,0,n) = -fine_state(i,j,0,n);
                        }

                        // Then, add the remaining positive proportionally.
                        if (SumP > 0.0) {
                            if (fine_state(i,j,0,n) > 0.0) {
                                Real alpha = (crseTot - std::abs(SumN)) / SumP;
                                fine(i,j,0,n) = alpha * fine_state(i,j,0,n);
                            }
                        } else { /* (SumP < 0) */
                            Real posVal = (crseTot - std::abs(SumN)) / cvol;
                            fine(i,j,0,n) += posVal;
                        }

                    }
                }

            } else if ( (crseTot > 0) && (crseTot < std::abs(SumN)) ) {

                /*
                 * Special case 2:
                 *
                 * Coarse correction > 0, and correction can not make
                 * them all positive.
                 *
                 * Add correction only to the negative cells
                 * in proportion to their magnitude, and
                 * don't add any correction to the states already positive.
                 */
                for     (int j = jlo; j <= jhi; ++j) {
                    for (int i = ilo; i <= ihi; ++i) {

                        Real alpha = crseTot / std::abs(SumN);
                        if (fine_state(i,j,0,n) < 0.0) {
                            // Add correction to negative cells proportionally.
                            fine(i,j,0,n) = alpha * std::abs(fine_state(i,j,0,n));
                        } else {
                            // Don't touch the positive states.
                            fine(i,j,0,n) = 0.0;
                        }

                    }
                }

            } else if ( (crseTot < 0) && (std::abs(crseTot) > SumP) ) {

                /*
                 * Special case 3:
                 *
                 * Coarse correction < 0, and fine_state DOES NOT have
                 * enough positive states to absorb it.
                 *
                 * Here we distribute the remaining negative amount
                 * in such a way as to make them all as close to the
                 * same negative value as possible.
                 */
                for     (int j = jlo; j <= jhi; ++j) {
                    for (int i = ilo; i <= ihi; ++i) {

                        // We want to make them all as close to the same negative value.
                        Real negVal = (SumP + SumN + crseTot) / cvol;
                        fine(i,j,0,n) = negVal - fine_state(i,j,0,n);

                    }
                }

            } else if ( (crseTot < 0) && (std::abs(crseTot) < SumP) &&
                        ((SumP+SumN+crseTot) > 0.0) )  {

                /*
                 * Special case 4:
                 *
                 * Coarse correction < 0, and fine_state has enough
                 * positive states to absorb all the negative
                 * correction *and* to redistribute to make
                 * negative cells positive.
                 */
                for     (int j = jlo; j <= jhi; ++j) {
                    for (int i = ilo; i <= ihi; ++i) {

                        if (fine_state(i,j,0,n) < 0.0) {
                            // Absorb the negative correction
                            fine(i,j,0,n) = -fine_state(i,j,0,n);
                        } else {
                            // Redistribute the rest to make negative cells positive
                            Real alpha = (crseTot + SumN) / SumP;
                            fine(i,j,0,n) = alpha * fine_state(i,j,0,n);
                        }

                    }
                }

            } else if ( (crseTot < 0) && (std::abs(crseTot) < SumP) &&
                        ((SumP+SumN+crseTot) < 0.0) )  {
                /*
                 * Special case 5:
                 *
                 * Coarse correction < 0, and fine_state has enough
                 * positive states to absorb all the negative
                 * correction, but not enough to fix the states
                 * already negative.
                 *
                 * Here we take a constant percentage away from each
                 * positive cell and don't touch the negatives.
                 */
                for     (int j = jlo; j <= jhi; ++j) {
                    for (int i = ilo; i <= ihi; ++i) {

                        if (fine_state(i,j,0,n) > 0.0) {
                            // Don't touch the negatives
                            fine(i,j,0,n) = -fine_state(i,j,0,n);
                        } else {
                            // Take a constant percentage away from each positive cell
                            Real alpha = (crseTot + SumP) / SumN;
                            fine(i,j,0,n) = alpha * fine_state(i,j,0,n);
                        }

                    }
                }

            } // special cases
        } // redo_me
    } // n

    // Set sync for density (n=0) to sum of spec sync (1:nvar)
    for     (int j = jlo; j <= jhi; ++j) {
        for (int i = ilo; i <= ihi; ++i) {
            fine(i,j,0,0) = 0.0;
            for (int n = 1; n < nvar-1; ++n) {
                fine(i,j,0,0) += fine(i,j,0,n);
            }
        }
    }

}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void ccquartic_interp (int i, int j, int /*k*/, int n,
                       Array4<Real const> const& crse,
                       Array4<Real>       const& fine) noexcept
{
    // Note: there are asserts in CellConservativeQuartic::interp()
    //       to check whether ratio is all equal to 2.

    constexpr Array1D<Real, -2, 2> cL = { -0.01171875_rt,  0.0859375_rt, 0.5_rt, -0.0859375_rt, 0.01171875_rt };

    int ic = amrex::coarsen(i,2);
    int jc = amrex::coarsen(j,2);
    int irx = i - 2*ic; // = abs(i % 2)
    int jry = j - 2*jc; // = abs(j % 2);

    Array1D<Real, -2, 2> ctmp;
    for (int ii = -2; ii <= 2; ++ii) {
        ctmp(ii) = 2.0_rt * ( cL(-2)*crse(ic+ii,jc-2,0,n)
                            + cL(-1)*crse(ic+ii,jc-1,0,n)
                            + cL( 0)*crse(ic+ii,jc,  0,n)
                            + cL( 1)*crse(ic+ii,jc+1,0,n)
                            + cL( 2)*crse(ic+ii,jc+2,0,n) );
        if (jry) {
            ctmp(ii) = 2.0_rt * crse(ic+ii,jc,0,n) - ctmp(ii);
        }
    } // ii

    Real ftmp = 2.0_rt * ( cL(-2)*ctmp(-2)
                         + cL(-1)*ctmp(-1)
                         + cL( 0)*ctmp( 0)
                         + cL( 1)*ctmp( 1)
                         + cL( 2)*ctmp( 2) );
    if (irx) {
        ftmp = 2.0_rt * ctmp(0) - ftmp;
    }

    fine(i,j,0,n) = ftmp;

}

}  // namespace amrex

#endif
