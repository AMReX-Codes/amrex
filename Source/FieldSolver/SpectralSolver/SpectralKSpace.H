#ifndef WARPX_SPECTRAL_K_SPACE_H_
#define WARPX_SPECTRAL_K_SPACE_H_

#include <AMReX_BoxArray.H>
#include <AMReX_LayoutData.H>

using SpectralKVector = amrex::LayoutData< amrex::Gpu::ManagedVector <amrex::Real> >;

/* TODO Documentation: class represent k space, and how it is distribution
 * for local FFT on each MPI: k spaces are not connected.
 */
class SpectralKSpace
{
    public:
        SpectralKSpace( const amrex::BoxArray& realspace_ba,
                        const amrex::DistributionMapping& dm,
                        const amrex::RealVect dx );
        SpectralKVector& AllocateAndFillKVector(
            const amrex::DistributionMapping& dm, const int i_dim ) const;
        SpectralKVector& AllocateAndFillModifiedKVector(
            const amrex::DistributionMapping& dm, const int i_dim, const int order ) const;
        amrex::BoxArray spectralspace_ba;

    private:
        amrex::Array<SpectralKVector, AMREX_SPACEDIM> k_vec;
        // 3D: array of 3 components, corresponding to kx, ky, kz
        // 2D: array of 2 components, corresponding to kx, kz
        amrex::RealVect dx;
};

#endif
