/* Copyright 2019 Axel Huebl
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef INTERPOLATE_DENSITY_FINE_TO_COARSE_H
#define INTERPOLATE_DENSITY_FINE_TO_COARSE_H

#include <AMReX_Array4.H>
#include <AMReX_REAL.H>
#include <AMReX_BLassert.H>
#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>

#include <utility> // std::move


/** Fill a charge density (rho) coarse patch with averaged values from a fine patch
 *
 * Fills the values of the charge density on the coarse patch
 * by averaging the values of the charge density of the fine patch.
 */
class InterpolateDensityFineToCoarse
{
public:
    /** Construct with fine and coarse patch and their refinement ratio
     *
     * @param[in]     fine read-only fine patch
     * @param[in,out] coarse overwritten coarse patch
     * @param[in]     refinement_ratio ratio between coarse and fine patch granularity
     *                (currently, only a value of is implemented)
     * @param[in]     number_of_components the number of components
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    InterpolateDensityFineToCoarse(
        amrex::Array4<amrex::Real const> const fine,
        amrex::Array4<amrex::Real      > const coarse,
        int const refinement_ratio,
        int const number_of_components
    ) : m_fine(std::move(fine)),
        m_coarse(std::move(coarse)),
        m_refinement_ratio(std::move(refinement_ratio)),
        m_number_of_components(std::move(number_of_components))
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
        auto const fine = [fine_unsafe] AMREX_GPU_DEVICE (int const j, int const k, int const l, int const m) noexcept
        {
            return fine_unsafe.contains(j, k, l) ? fine_unsafe(j, k, l, m) :  amrex::Real{0.};
        };

        int const ii = i * m_refinement_ratio;
        int const jj = j * m_refinement_ratio;
        int const kk = k * m_refinement_ratio;
        for( int m = 0; m < m_number_of_components; ++m ) {
#if AMREX_SPACEDIM == 2
            coarse(i,j,kk,m) = 0.25_rt * (
                fine(ii,jj,kk,m) +
                0.5_rt  * (
                    fine(ii-1,jj  ,kk,m) + fine(ii+1,jj  ,kk,m) +
                    fine(ii  ,jj-1,kk,m) + fine(ii  ,jj+1,kk,m)
                ) +
                0.25_rt * (
                    fine(ii-1,jj-1,kk,m) + fine(ii+1,jj-1,kk,m) +
                    fine(ii-1,jj+1,kk,m) + fine(ii+1,jj+1,kk,m)
                )
            );
#elif AMREX_SPACEDIM == 3
            coarse(i,j,k,m) = 0.125_rt * (
                fine(ii,jj,kk,m) +
                0.5_rt   * (
                    fine(ii-1,jj  ,kk  ,m) + fine(ii+1,jj  ,kk  ,m) +
                    fine(ii  ,jj-1,kk  ,m) + fine(ii  ,jj+1,kk  ,m) +
                    fine(ii  ,jj  ,kk-1,m) + fine(ii  ,jj  ,kk+1,m)
                ) +
                0.25_rt  * (
                    fine(ii-1,jj-1,kk  ,m) + fine(ii+1,jj-1,kk  ,m) +
                    fine(ii-1,jj+1,kk  ,m) + fine(ii+1,jj+1,kk  ,m) +
                    fine(ii-1,jj  ,kk-1,m) + fine(ii+1,jj  ,kk-1,m) +
                    fine(ii-1,jj  ,kk+1,m) + fine(ii+1,jj  ,kk+1,m) +
                    fine(ii  ,jj-1,kk-1,m) + fine(ii  ,jj+1,kk-1,m) +
                    fine(ii  ,jj-1,kk+1,m) + fine(ii  ,jj+1,kk+1,m)
                ) +
                0.125_rt * (
                    fine(ii-1,jj-1,kk-1,m) + fine(ii-1,jj-1,kk+1,m) +
                    fine(ii-1,jj+1,kk-1,m) + fine(ii-1,jj+1,kk+1,m) +
                    fine(ii+1,jj-1,kk-1,m) + fine(ii+1,jj-1,kk+1,m) +
                    fine(ii+1,jj+1,kk-1,m) + fine(ii+1,jj+1,kk+1,m)
                )
            );
#endif
        }
    }
private:
    amrex::Array4< amrex::Real const > m_fine;
    amrex::Array4< amrex::Real       > m_coarse;
    int m_refinement_ratio;
    int m_number_of_components;
};

#endif // INTERPOLATE_DENSITY_FINE_TO_COARSE_H
