#ifndef WARPX_SPECTRAL_FIELD_DATA_H_
#define WARPX_SPECTRAL_FIELD_DATA_H_

#include <WarpX_ComplexForFFT.H>
#include <SpectralKSpace.H>
#include <AMReX_MultiFab.H>

// Declare type for spectral fields
using SpectralField = amrex::FabArray< amrex::BaseFab <Complex> >;

/* Index for the regular fields, when stored in spectral space */
struct SpectralFieldIndex {
  enum { Ex=0, Ey, Ez, Bx, By, Bz, Jx, Jy, Jz, rho_old, rho_new, n_fields };
  // n_fields is automatically the total number of fields
};

/* Index for the PML fields, when stored in spectral space */
struct SpectralPMLIndex {
  enum { Exy=0, Exz, Eyx, Eyz, Ezx, Ezy,
           Bxy, Bxz, Byx, Byz, Bzx, Bzy, n_fields };
  // n_fields is automatically the total number of fields
};

/* \brief Class that stores the fields in spectral space, and performs the
 *  Fourier transforms between real space and spectral space
 */
class SpectralFieldData
{

    // Define the FFTplans type, which holds one fft plan per box
    // (plans are only initialized for the boxes that are owned by
    // the local MPI rank)
#ifdef AMREX_USE_GPU
    using FFTplans = amrex::LayoutData<cufftHandle>;
#else
    using FFTplans = amrex::LayoutData<fftw_plan>;
#endif

    public:
        SpectralFieldData( const amrex::BoxArray& realspace_ba,
                           const SpectralKSpace& k_space,
                           const amrex::DistributionMapping& dm,
                           const int n_field_required );
        SpectralFieldData() = default; // Default constructor
        SpectralFieldData& operator=(SpectralFieldData&& field_data) = default;
        ~SpectralFieldData();
        void ForwardTransform( const amrex::MultiFab& mf,
                               const int field_index, const int i_comp);
        void BackwardTransform( amrex::MultiFab& mf,
                               const int field_index, const int i_comp);
        // `fields` stores fields in spectral space, as multicomponent FabArray
        SpectralField fields;

    private:
        // tmpRealField and tmpSpectralField store fields
        // right before/after the Fourier transform
        SpectralField tmpSpectralField; // contains Complexs
        amrex::MultiFab tmpRealField; // contains Reals
        FFTplans forward_plan, backward_plan;
        // Correcting "shift" factors when performing FFT from/to
        // a cell-centered grid in real space, instead of a nodal grid
        SpectralShiftFactor xshift_FFTfromCell, xshift_FFTtoCell,
                            zshift_FFTfromCell, zshift_FFTtoCell;
#if (AMREX_SPACEDIM==3)
        SpectralShiftFactor yshift_FFTfromCell, yshift_FFTtoCell;
#endif
};

#endif // WARPX_SPECTRAL_FIELD_DATA_H_
