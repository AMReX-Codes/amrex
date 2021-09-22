#ifndef AMREX_INTERP_C_H_
#define AMREX_INTERP_C_H_
#include <AMReX_Config.H>
#include <AMReX_Array4.H>

#if (AMREX_SPACEDIM == 1)
#include <AMReX_Interp_1D_C.H>
#elif (AMREX_SPACEDIM == 2)
#include <AMReX_Interp_2D_C.H>
#else
#include <AMReX_Interp_3D_C.H>
#endif

namespace amrex {

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void ccprotect_redo (int ic, int jc, int kc, int n,
                     Box const& fine_bx,
                     IntVect const& ratio,
#if (AMREX_SPACEDIM == 2)
                     GeometryData cs_geomdata,
                     GeometryData fn_geomdata,
#endif
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
#if (AMREX_SPACEDIM == 2)
    int klo = 0;
    int khi = 0;
#else
    int klo = amrex::max(ratio[2]*kc,              fnbxlo.z);
    int khi = amrex::min(ratio[2]*kc+(ratio[2]-1), fnbxhi.z);
#endif

#if (AMREX_SPACEDIM == 2)
    /*
     * Calculate coarse cell volumes.
     */
    Real cvol;
    ccprotect_calc_cvol(cvol, ic, cs_geomdata);
#else
    // Calculate number of fine cells
    int numFineCells = (ihi-ilo+1) * (jhi-jlo+1) * (khi-klo+1);
#endif

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
    for         (int k = klo; k <= khi; ++k) {
        for     (int j = jlo; j <= jhi; ++j) {
            for (int i = ilo; i <= ihi; ++i) {
                ccprotect_calc_sums(i, j, k, n,
                                    crseTot, SumN, SumP,
#if (AMREX_SPACEDIM == 2)
                                    fn_geomdata,
#endif
                                    fine, fine_state);
            }
        }
    }

    if ( (crseTot > 0) && (crseTot > Math::abs(SumN)) ) {

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
        for         (int k = klo; k <= khi; ++k) {
            for     (int j = jlo; j <= jhi; ++j) {
                for (int i = ilo; i <= ihi; ++i) {
                    ccprotect_case1(i, j, k, n,
                                    crseTot, SumN, SumP,
#if (AMREX_SPACEDIM == 2)
                                    cvol,
#else
                                    numFineCells,
#endif
                                    fine, fine_state);
                }
            }
        }

    } else if ( (crseTot > 0) && (crseTot < Math::abs(SumN)) ) {

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
        for         (int k = klo; k <= khi; ++k) {
            for     (int j = jlo; j <= jhi; ++j) {
                for (int i = ilo; i <= ihi; ++i) {
                    ccprotect_case2(i, j, k, n,
                                    crseTot, SumN, SumP,
                                    fine, fine_state);
                }
            }
        }

    } else if ( (crseTot < 0) && (Math::abs(crseTot) > SumP) ) {

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
        for         (int k = klo; k <= khi; ++k) {
            for     (int j = jlo; j <= jhi; ++j) {
                for (int i = ilo; i <= ihi; ++i) {
                    ccprotect_case3(i, j, k, n,
                                    crseTot, SumN, SumP,
#if (AMREX_SPACEDIM == 2)
                                    cvol,
#else
                                    numFineCells,
#endif
                                    fine, fine_state);
                }
            }
        }

    } else if ( (crseTot < 0) && (Math::abs(crseTot) < SumP) &&
                ((SumP+SumN+crseTot) > 0.0) )  {

        /*
         * Special case 4:
         *
         * Coarse correction < 0, and fine_state has enough
         * positive states to absorb all the negative
         * correction *and* to redistribute to make
         * negative cells positive.
         */
        for         (int k = klo; k <= khi; ++k) {
            for     (int j = jlo; j <= jhi; ++j) {
                for (int i = ilo; i <= ihi; ++i) {
                    ccprotect_case4(i, j, k, n,
                                    crseTot, SumN, SumP,
                                    fine, fine_state);
                }
            }
        }

    } else if ( (crseTot < 0) && (Math::abs(crseTot) < SumP) &&
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
        for         (int k = klo; k <= khi; ++k) {
            for     (int j = jlo; j <= jhi; ++j) {
                for (int i = ilo; i <= ihi; ++i) {
                    ccprotect_case5(i, j, k, n,
                                    crseTot, SumN, SumP,
                                    fine, fine_state);
                }
            }
        }

    }

}

}  // namespace amrex

#endif
