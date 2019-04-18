#ifndef WARPX_PSATD_DATA_H_
#define WARPX_PSATD_DATA_H_

#include <WarpX_Complex.H>
#include <AMReX_MultiFab.H>

using namespace amrex;

// Declare type for spectral fields
using SpectralField = FabArray<BaseFab<Complex>>;

/* Fields that will be stored in spectral space */
struct SpectralFieldIndex{
    enum { Ex=0, Ey, Ez, Bx, By, Bz, Jx, Jy, Jz, rho_old, rho_new };
};

/* \brief Class that stores the fields in spectral space
 *  and performs the spectral transforms to/from real space
 */
class SpectralData
{
    friend class PsatdSolver;

#ifdef AMREX_USE_GPU
// Add cuFFT-specific code
#else
    using FFTplans = LayoutData<fftw_plan>;
#endif

    public:
        SpectralData( const BoxArray& realspace_ba,
                      const BoxArray& spectralspace_ba,
                      const DistributionMapping& dm );
        ~SpectralData();
        void ForwardTransform( const MultiFab& mf, const int field_index );
        void BackwardTransform( MultiFab& mf, const int field_index );

    private:
        SpectralField Ex, Ey, Ez, Bx, By, Bz, Jx, Jy, Jz, rho_old, rho_new;
        SpectralField tmpRealField, tmpSpectralField; // Store fields before/after transform
        FFTplans forward_plan, backward_plan;
        SpectralField& getSpectralField( const int field_index );
};

#endif // WARPX_PSATD_DATA_H_
