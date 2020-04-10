#ifndef WARPX_AVERAGE_H_
#define WARPX_AVERAGE_H_

#include "WarpX.H"

namespace Average{

    using namespace amrex;

    /**
     * \brief Interpolates the floating point data contained in the source Array4
     *        \c mf_fp_arr, extracted from a fine MultiFab, by averaging over either
     *        1 point or 2 equally distant points.
     *
     * \param[in] arr_src floating point data to be interpolated
     * \param[in] sf      staggering of the source fine MultiFab
     * \param[in] sc      staggering of the destination coarsened MultiFab
     * \param[in] cr      coarsening ratio along each spatial direction
     * \param[in] np      number of points to loop over for interpolation
     * \param[in] i       index along x of the coarsened Array4 to be filled
     * \param[in] j       index along y of the coarsened Array4 to be filled
     * \param[in] k       index along z of the coarsened Array4 to be filled
     * \param[in] comp    index along the fourth component of the Array4 \c arr_src
     *                    containing the data to be interpolated
     *
     * \return interpolated field at cell (i,j,k) of a coarsened Array4
     */
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    Real CoarsenAndInterpolateKernel ( Array4<Real const> const& arr_src,
                                       int const* const AMREX_RESTRICT sf,
                                       int const* const AMREX_RESTRICT sc,
                                       int const* const AMREX_RESTRICT cr,
                                       int const* const AMREX_RESTRICT np,
                                       const int i,
                                       const int j,
                                       const int k,
                                       const int comp )
    {
        // Indices of source coarse array
        const int ic[3] = { i, j, k };

        // Index of first point of source Array4 from which interpolation is done
        int idx_min[3];

        // Compute number of points to loop over (either 1 or 2)
        // and starting indices of fine array in each direction
        for ( int l = 0; l < 3; ++l ) {
            if ( cr[l] == 1 ) idx_min[l] = ic[l]-sc[l]*(1-sf[l]); // no coarsening
            else              idx_min[l] = ic[l]*cr[l]+static_cast<int>(cr[l]/2)*(1-sc[l])-(1-sf[l]);
        }

        // Auxiliary integer variables
        const int numx = np[0];
        const int numy = np[1];
        const int numz = np[2];
        const int imin = idx_min[0];
        const int jmin = idx_min[1];
        const int kmin = idx_min[2];

        // Interpolate over points computed above
        Real c = 0.0_rt;
        for (int kref = 0; kref < numz; ++kref) {
            for (int jref = 0; jref < numy; ++jref) {
                for (int iref = 0; iref < numx; ++iref) {
                    c += arr_src(imin+iref,jmin+jref,kmin+kref,comp);
                }
            }
        }
        return c / (numx*numy*numz);
    };

    /**
     * \brief Loops over the boxes of the coarsened MultiFab \c mf_dst and fills
     *        them by interpolating the data contained in the fine MultiFab \c mf_src.
     *
     * \param[in,out] mf_dst     coarsened MultiFab containing the floating point data
     *                           to be filled by interpolating the source fine MultiFab
     * \param[in]     mf_src     fine MultiFab containing the floating point data to be interpolated
     * \param[in]     dcomp      offset for the fourth component of the coarsened Array4
     *                           object, extracted from its MultiFab, where the interpolated
     *                           values will be stored
     * \param[in]     scomp      offset for the fourth component of the fine Array4
     *                           object, extracted from its MultiFab, containing the
     *                           data to be interpolated
     * \param[in]     ncomp      number of components to loop over for the coarsened
     *                           Array4 extracted from the coarsened MultiFab \c mf_dst
     * \param[in]     ngrow      number of guard cells to fill
     * \param[in]     crse_ratio coarsening ratio between the fine MultiFab \c mf_src
     *                           and the coarsened MultiFab \c mf_dst
     */
    void CoarsenAndInterpolateLoop ( MultiFab& mf_dst,
                                     const MultiFab& mf_src,
                                     const int dcomp,
                                     const int scomp,
                                     const int ncomp,
                                     const int ngrow,
                                     const IntVect crse_ratio=IntVect(1) );

    /**
     * \brief Stores in the coarsened MultiFab \c mf_dst the values obtained by
     *        interpolating the data contained in the fine MultiFab \c mf_src.
     *
     * \param[in,out] mf_dst     coarsened MultiFab containing the floating point data
     *                           to be filled by interpolating the fine MultiFab \c mf_src
     * \param[in]     mf_src     fine MultiFab containing the floating point data to be interpolated
     * \param[in]     dcomp      offset for the fourth component of the coarsened Array4
     *                           object, extracted from its MultiFab, where the interpolated
     *                           values will be stored
     * \param[in]     scomp      offset for the fourth component of the fine Array4
     *                           object, extracted from its MultiFab, containing the
     *                           data to be interpolated
     * \param[in]     ncomp      number of components to loop over for the coarsened
     *                           Array4 extracted from the coarsened MultiFab \c mf_dst
     * \param[in]     ngrow      number of guard cells to fill
     * \param[in]     crse_ratio coarsening ratio between the fine MultiFab \c mf_src
     *                           and the coarsened MultiFab \c mf_dst
     */
    void CoarsenAndInterpolate ( MultiFab& mf_dst,
                                 const MultiFab& mf_src,
                                 const int dcomp,
                                 const int scomp,
                                 const int ncomp,
                                 const int ngrow,
                                 const IntVect crse_ratio=IntVect(1) );
}

#endif // WARPX_AVERAGE_H_
