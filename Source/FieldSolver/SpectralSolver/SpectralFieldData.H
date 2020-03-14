/* Copyright 2019 David Grote, Maxence Thevenet, Remi Lehe
 * Revathi Jambunathan
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_SPECTRAL_FIELD_DATA_H_
#define WARPX_SPECTRAL_FIELD_DATA_H_

#include "Utils/WarpX_Complex.H"
#include "SpectralKSpace.H"
#include <AMReX_MultiFab.H>

#include <string>

// Declare type for spectral fields
using SpectralField = amrex::FabArray< amrex::BaseFab <Complex> >;

/** Index for the regular fields, when stored in spectral space:
 *  - n_fields is automatically the total number of fields
 *  - divE reuses the memory slot for Bx, since Bx is not used when computing divE
 */
struct SpectralFieldIndex {
  enum { Ex=0, Ey, Ez, Bx, By, Bz, Jx, Jy, Jz, rho_old, rho_new, n_fields, divE=3 };
};

/** Index for the PML fields, when stored in spectral space */
struct SpectralPMLIndex {
  enum { Exy=0, Exz, Eyx, Eyz, Ezx, Ezy,
           Bxy, Bxz, Byx, Byz, Bzx, Bzy, n_fields };
  // n_fields is automatically the total number of fields
};

/** \brief Class that stores the fields in spectral space, and performs the
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
#  ifdef AMREX_USE_FLOAT
    using FFTplans = amrex::LayoutData<fftwf_plan>;
#  else
    using FFTplans = amrex::LayoutData<fftw_plan>;
#  endif
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

#ifdef AMREX_USE_GPU
        /** \brief This method converts a cufftResult
        * into the corresponding string
        *
        * @param[in] err a cufftResult
        * @return an std::string
        */
        std::string cufftErrorToString (const cufftResult& err);
#endif
};

#endif // WARPX_SPECTRAL_FIELD_DATA_H_
