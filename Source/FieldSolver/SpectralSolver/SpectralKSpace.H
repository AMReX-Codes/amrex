#ifndef WARPX_SPECTRAL_K_SPACE_H_
#define WARPX_SPECTRAL_K_SPACE_H_

#include <WarpX_Complex.H>
#include <AMReX_BoxArray.H>
#include <AMReX_LayoutData.H>

using KVectorComponent = amrex::LayoutData< amrex::Gpu::ManagedVector <amrex::Real> >;
using SpectralShiftFactor = amrex::LayoutData< amrex::Gpu::ManagedVector <Complex> >;

struct ShiftType {
    enum{ NodalToCentered=0, CenteredToNodal=1 };
};

/* TODO Documentation: class represent k space, and how it is distribution
 * for local FFT on each MPI: k spaces are not connected.
 */
class SpectralKSpace
{
    public:
        amrex::BoxArray spectralspace_ba;
        SpectralKSpace( const amrex::BoxArray& realspace_ba,
                        const amrex::DistributionMapping& dm,
                        const amrex::RealVect dx );
        KVectorComponent AllocateAndFillKComponent(
            const amrex::DistributionMapping& dm, const int i_dim ) const;
        KVectorComponent AllocateAndFillModifiedKComponent(
            const amrex::DistributionMapping& dm, const int i_dim, const int order ) const;
        SpectralShiftFactor AllocateAndFillSpectralShiftFactor(
            const amrex::DistributionMapping& dm, const int i_dim, const int order ) const;


    private:
        amrex::Array<KVectorComponent, AMREX_SPACEDIM> k_vec;
        // 3D: array of 3 components, corresponding to kx, ky, kz
        // 2D: array of 2 components, corresponding to kx, kz
        amrex::RealVect dx;
};

#endif
