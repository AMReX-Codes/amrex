#ifndef WARPX_SPECTRAL_FIELD_DATA_H_
#define WARPX_SPECTRAL_FIELD_DATA_H_

#include <WarpX_Complex.H>
#include <SpectralKSpace.H>
#include <AMReX_MultiFab.H>

// Declare type for spectral fields
using SpectralField = amrex::FabArray< amrex::BaseFab <Complex> >;

/* Fields that will be stored in spectral space */
struct SpectralFieldIndex{
    enum { Ex=0, Ey, Ez, Bx, By, Bz, Jx, Jy, Jz, rho_old, rho_new };
};

/* \brief Class that stores the fields in spectral space
 *  and performs the spectral transforms to/from real space
 */
class SpectralFieldData
{
    friend class PsatdAlgorithm;

#ifdef AMREX_USE_GPU
// Add cuFFT-specific code
#else
    using FFTplans = amrex::LayoutData<fftw_plan>;
#endif

    public:
        SpectralFieldData( const amrex::BoxArray& realspace_ba,
                      const SpectralKSpace& k_space,
                      const amrex::DistributionMapping& dm );
        SpectralFieldData() = default; // Default constructor
        SpectralFieldData& operator=(SpectralFieldData&& field_data) = default; // Default move assignment
        ~SpectralFieldData();
        void ForwardTransform( const amrex::MultiFab& mf, const int field_index );
        void BackwardTransform( amrex::MultiFab& mf, const int field_index );

    private:
        SpectralField Ex, Ey, Ez, Bx, By, Bz, Jx, Jy, Jz, rho_old, rho_new;
        SpectralField tmpRealField, tmpSpectralField; // Store fields before/after transform
        FFTplans forward_plan, backward_plan;
        SpectralField& getSpectralField( const int field_index );
};

#endif // WARPX_SPECTRAL_FIELD_DATA_H_
