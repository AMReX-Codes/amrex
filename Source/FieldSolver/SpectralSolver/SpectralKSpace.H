#ifndef WARPX_SPECTRAL_K_SPACE_H_
#define WARPX_SPECTRAL_K_SPACE_H_

#include <AMReX_BoxArray.H>
#include <AMReX_LayoutData.H>

using namespace amrex;
using namespace Gpu;

using SpectralKVector = LayoutData<ManagedVector<Real>>;

/* TODO Documentation: class represent k space, and how it is distribution
 * for local FFT on each MPI: k spaces are not connected.
 */
class SpectralKSpace
{
    public:
        SpectralKSpace( const BoxArray& realspace_ba, const DistributionMapping& dm, const Real* dx );
        SpectralKVector& getModifiedKVector( const int i_dim ) const;
        BoxArray spectralspace_ba;

    private:
        SpectralKVector kx_vec, ky_vec, kz_vec;
        const Real* dx;
};

void
AllocateAndFillKvector( ManagedVector<Real>& k, const Box& bx, const Real* dx, const int i_dim );

#endif
