/* Copyright 2019-2020 Axel Huebl
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef INTERPOLATECURRENTFINETOCOARSE_H
#define INTERPOLATECURRENTFINETOCOARSE_H

#include <AMReX_Array4.H>
#include <AMReX_REAL.H>
#include <AMReX_BLassert.H>
#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>

#include <utility>  // std::move


/** Fill a current coarse patch with averaged values from a fine patch
 *
 * Fills the values of the current for a selected component on the coarse patch
 * by averaging the values of the current of the fine patch.
 *
 * @tparam IDim j-field component on which the averaging is performed
 */
template< int IDim >
class InterpolateCurrentFineToCoarse
{
public:
    /** Construct with fine and coarse patch and their refinement ratio
     *
     * @param[in]     fine read-only fine patch
     * @param[in,out] coarse overwritten coarse patch
     * @param[in]     refinement_ratio ratio between coarse and fine patch granularity
     *                (currently, only a value of is implemented)
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    InterpolateCurrentFineToCoarse(
        amrex::Array4< amrex::Real const > const fine,
        amrex::Array4< amrex::Real       > const coarse,
        int const refinement_ratio
    ) : m_fine(std::move(fine)),
        m_coarse(std::move(coarse)),
        m_refinement_ratio(std::move(refinement_ratio))
    {
        //! @note constants and stencils in operator() implementation assume 2x refinement
        BL_ASSERT(refinement_ratio == 2);
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void
    operator()(
        int const i,
        int const j,
        int const k
    ) const noexcept
    {
        using namespace amrex;

        auto const & fine_unsafe = m_fine;   // out-of-bounds access not secured with zero-values yet
        auto const & coarse      = m_coarse; // out-of-bounds access not secured but will also not occur

        // return zero for out-of-bounds accesses during interpolation
        // this is efficiently used as a method to add neutral elements beyond guards in the average below
        auto const fine = [fine_unsafe] AMREX_GPU_DEVICE (int const jj, int const kk, int const ll) noexcept
        {
            return fine_unsafe.contains(jj, kk, ll) ? fine_unsafe(jj, kk, ll) :  amrex::Real{0.};
        };

        int const ii = i * m_refinement_ratio;
        int const jj = j * m_refinement_ratio;
        int const kk = k * m_refinement_ratio;
#if AMREX_SPACEDIM == 2
        if (IDim == 0) {
            coarse(i, j, k) = 0.25_rt * (
                fine(ii, jj, kk) + fine(ii + 1, jj, kk) +
                0.5_rt  * (
                    fine(ii, jj - 1, kk) + fine(ii + 1, jj - 1, kk) +
                    fine(ii, jj + 1, kk) + fine(ii + 1, jj + 1, kk)
                )
            );
        } else if (IDim == 2) {
            coarse(i, j, k) = 0.25_rt * (
                fine(ii, jj, kk) + fine(ii, jj + 1, kk) +
                0.5_rt  * (
                    fine(ii - 1, jj, kk) + fine(ii - 1, jj + 1, kk) +
                    fine(ii + 1, jj, kk) + fine(ii + 1, jj + 1, kk)
                )
            );
        } else {
            coarse(i, j, k) = 0.25_rt * (
                fine(ii, jj, kk) +
                0.5_rt  * (
                    fine(ii - 1, jj    , kk) + fine(ii + 1, jj    , kk) +
                    fine(ii    , jj - 1, kk) + fine(ii    , jj + 1, kk)
                ) +
                0.25_rt * (
                    fine(ii - 1, jj - 1, kk) + fine(ii + 1, jj - 1, kk) +
                    fine(ii - 1, jj + 1, kk) + fine(ii + 1, jj + 1, kk)
                )
            );
        }
#elif AMREX_SPACEDIM == 3
        if (IDim == 0) {
            coarse(i,j,k) = 0.125_rt * (
                fine(ii  , jj, kk) +
                0.5_rt  * (
                    fine(ii  , jj-1, kk  ) + fine(ii  , jj+1, kk  ) +
                    fine(ii  , jj  , kk-1) + fine(ii  , jj  , kk+1)
                ) +
                0.25_rt * (
                    fine(ii  , jj-1, kk-1) + fine(ii  , jj+1, kk-1) +
                    fine(ii  , jj-1, kk+1) + fine(ii  , jj+1, kk+1)
                ) +
                fine(ii+1, jj, kk) +
                0.5_rt  * (
                    fine(ii+1, jj-1, kk  ) + fine(ii+1, jj+1, kk  ) +
                    fine(ii+1, jj  , kk-1) + fine(ii+1, jj  , kk+1)
                ) +
                0.25_rt * (
                    fine(ii+1, jj-1, kk-1) + fine(ii+1, jj+1, kk-1) +
                    fine(ii+1, jj-1, kk+1) + fine(ii+1, jj+1, kk+1)
                )
            );
        } else if (IDim == 1) {
            coarse(i, j, k) = 0.125_rt * (
                fine(ii, jj , kk) +
                0.5_rt  * (
                    fine(ii-1, jj  , kk  ) + fine(ii+1, jj  , kk  ) +
                    fine(ii  , jj  , kk-1) + fine(ii  , jj  , kk+1)
                ) +
                0.25_rt * (
                    fine(ii-1, jj  , kk-1) + fine(ii+1, jj  , kk-1) +
                    fine(ii-1, jj  , kk+1) + fine(ii+1, jj  , kk+1)
                ) +
                fine(ii, jj+1, kk) +
                0.5_rt  * (
                    fine(ii-1, jj+1, kk  ) + fine(ii+1, jj+1, kk  ) +
                    fine(ii  , jj+1, kk-1) + fine(ii  , jj+1, kk+1)
                ) +
                0.25_rt * (
                    fine(ii-1, jj+1, kk-1) + fine(ii+1, jj+1, kk-1) +
                    fine(ii-1, jj+1, kk+1) + fine(ii+1, jj+1, kk+1)
                )
             );
        } else {
            coarse(i, j, k) = 0.125_rt * (
                fine(ii, jj, kk  ) +
                0.5_rt  * (
                    fine(ii-1, jj  , kk  ) + fine(ii+1, jj  , kk  ) +
                    fine(ii  , jj-1, kk  ) + fine(ii  , jj+1, kk  )
                ) +
                0.25_rt * (
                    fine(ii-1, jj-1, kk  ) + fine(ii+1, jj-1, kk  ) +
                    fine(ii-1, jj+1, kk  ) + fine(ii+1, jj+1, kk  )
                ) +
                fine(ii, jj, kk+1) +
                0.5_rt  * (
                    fine(ii-1, jj  , kk+1) + fine(ii+1, jj  , kk+1) +
                    fine(ii  , jj-1, kk+1) + fine(ii  , jj+1, kk+1)
                ) +
                0.25_rt * (
                    fine(ii-1, jj-1, kk+1) + fine(ii+1, jj-1, kk+1) +
                    fine(ii-1, jj+1, kk+1) + fine(ii+1, jj+1, kk+1)
                )
             );
         }
#endif
    }
private:
    amrex::Array4< amrex::Real const > m_fine;
    amrex::Array4< amrex::Real       > m_coarse;
    int m_refinement_ratio;
};

#endif // INTERPOLATECURRENTFINETOCOARSE_H
