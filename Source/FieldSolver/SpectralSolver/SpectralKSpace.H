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
        amrex::BoxArray spectralspace_ba;
        SpectralKVector kx_vec, kz_vec;
#if (AMREX_SPACEDIM==3)
        SpectralKVector ky_vec;
#endif
        amrex::RealVect dx;
};

void
AllocateAndFillKvector( amrex::Gpu::ManagedVector<amrex::Real>& k,
    const amrex::Box& bx, const amrex::RealVect dx, const int i_dim );

void
ComputeModifiedKVector( amrex::Gpu::ManagedVector<amrex::Real>& modified_k,
                        const amrex::Gpu::ManagedVector<amrex::Real>& k,
                        const amrex::Box& bx, const amrex::Real dx,
                        const int norder );

#endif
