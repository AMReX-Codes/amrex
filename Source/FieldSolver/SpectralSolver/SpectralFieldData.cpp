#include <SpectralFieldData.H>

using namespace amrex;

SpectralFieldData::SpectralFieldData( const BoxArray& realspace_ba,
                            const SpectralKSpace& k_space,
                            const DistributionMapping& dm )
{
    const BoxArray& spectralspace_ba = k_space.spectralspace_ba;

    // Allocate the arrays that contain the fields in spectral space
    Ex = SpectralField(spectralspace_ba, dm, 1, 0);
    Ey = SpectralField(spectralspace_ba, dm, 1, 0);
    Ez = SpectralField(spectralspace_ba, dm, 1, 0);
    Bx = SpectralField(spectralspace_ba, dm, 1, 0);
    By = SpectralField(spectralspace_ba, dm, 1, 0);
    Bz = SpectralField(spectralspace_ba, dm, 1, 0);
    Jx = SpectralField(spectralspace_ba, dm, 1, 0);
    Jy = SpectralField(spectralspace_ba, dm, 1, 0);
    Jz = SpectralField(spectralspace_ba, dm, 1, 0);
    rho_old = SpectralField(spectralspace_ba, dm, 1, 0);
    rho_new = SpectralField(spectralspace_ba, dm, 1, 0);

    // Allocate temporary arrays - over different boxarrays
    tmpRealField = SpectralField(realspace_ba, dm, 1, 0);
    tmpSpectralField = SpectralField(spectralspace_ba, dm, 1, 0);

    // Allocate and initialize the FFT plans
    forward_plan = FFTplans(spectralspace_ba, dm);
    backward_plan = FFTplans(spectralspace_ba, dm);
    for ( MFIter mfi(spectralspace_ba, dm); mfi.isValid(); ++mfi ){
        Box bx = spectralspace_ba[mfi];
#ifdef AMREX_USE_GPU
        // Add cuFFT-specific code
#else
        // Create FFTW plans
        forward_plan[mfi] = fftw_plan_dft_3d(
            // Swap dimensions: AMReX data is Fortran-order, but FFTW is C-order
            bx.length(2), bx.length(1), bx.length(0),
            reinterpret_cast<fftw_complex*>( tmpRealField[mfi].dataPtr() ),
            reinterpret_cast<fftw_complex*>( tmpSpectralField[mfi].dataPtr() ),
            FFTW_FORWARD, FFTW_ESTIMATE );
        backward_plan[mfi] = fftw_plan_dft_3d(
            // Swap dimensions: AMReX data is Fortran-order, but FFTW is C-order
            bx.length(2), bx.length(1), bx.length(0),
            reinterpret_cast<fftw_complex*>( tmpSpectralField[mfi].dataPtr() ),
            reinterpret_cast<fftw_complex*>( tmpRealField[mfi].dataPtr() ),
            FFTW_BACKWARD, FFTW_ESTIMATE );
        // TODO: Add 2D code
        // TODO: Do real-to-complex transform instead of complex-to-complex
        // TODO: Let the user decide whether to use FFTW_ESTIMATE or FFTW_MEASURE
#endif
    }
}


SpectralFieldData::~SpectralFieldData()
{
    if (tmpRealField.size() > 0){
        for ( MFIter mfi(tmpRealField); mfi.isValid(); ++mfi ){
#ifdef AMREX_USE_GPU
            // Add cuFFT-specific code
#else
            // Destroy FFTW plans
            fftw_destroy_plan( forward_plan[mfi] );
            fftw_destroy_plan( backward_plan[mfi] );
#endif
        }
    }
}

/* TODO: Documentation
 * Example: ForwardTransform( Efield_cp[0], SpectralFieldIndex::Ex )
 */
void
SpectralFieldData::ForwardTransform( const MultiFab& mf, const int field_index )
{
    // Loop over boxes
    for ( MFIter mfi(mf); mfi.isValid(); ++mfi ){

        // Copy the real-space field `mf` to the temporary field `tmpRealField`
        // This ensures that all fields have the same number of points
        // before the Fourier transform.
        // As a consequence, the copy discards the *last* point of `mf`
        // in any direction that has *nodal* index type.
        {
            Box bx = mf[mfi].box();
            const Box realspace_bx = bx.enclosedCells(); // discards last point in each nodal direction
            AMREX_ALWAYS_ASSERT( realspace_bx == tmpRealField[mfi].box() );
            Array4<const Real> mf_arr = mf[mfi].array();
            Array4<Complex> tmp_arr = tmpRealField[mfi].array();
            ParallelFor( realspace_bx,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                tmp_arr(i,j,k) = mf_arr(i,j,k);
            });
        }

        // Perform Fourier transform from `tmpRealField` to `tmpSpectralField`
#ifdef AMREX_USE_GPU
        // Add cuFFT-specific code ; make sure that this is done on the same
        // GPU stream as the above copy
#else
        fftw_execute( forward_plan[mfi] );
#endif

        // Copy the spectral-space field `tmpSpectralField` to the appropriate field
        // (specified by the input argument field_index )
        {
            SpectralField& field = getSpectralField( field_index );
            Array4<Complex> field_arr = field[mfi].array();
            Array4<const Complex> tmp_arr = tmpSpectralField[mfi].array();
            const Box spectralspace_bx = tmpSpectralField[mfi].box();
            ParallelFor( spectralspace_bx,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                field_arr(i,j,k) = tmp_arr(i,j,k);
            });
        }
    }
}


/* TODO: Documentation
 */
void
SpectralFieldData::BackwardTransform( MultiFab& mf, const int field_index )
{
    // Loop over boxes
    for ( MFIter mfi(mf); mfi.isValid(); ++mfi ){

        // Copy the appropriate field (specified by the input argument field_index)
        // to the spectral-space field `tmpSpectralField`
        {
            SpectralField& field = getSpectralField( field_index );
            Array4<const Complex> field_arr = field[mfi].array();
            Array4<Complex> tmp_arr = tmpSpectralField[mfi].array();
            const Box spectralspace_bx = tmpSpectralField[mfi].box();
            ParallelFor( spectralspace_bx,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                tmp_arr(i,j,k) = field_arr(i,j,k);
            });
        }

        // Perform Fourier transform from `tmpSpectralField` to `tmpRealField`
#ifdef AMREX_USE_GPU
        // Add cuFFT-specific code ; make sure that this is done on the same
        // GPU stream as the above copy
#else
        fftw_execute( backward_plan[mfi] );
#endif

        // Copy the temporary field `tmpRealField` to the real-space field `mf`
        // The copy does *not* fill the *last* point of `mf`
        // in any direction that has *nodal* index type (but this point is
        // in the guard cells and will be filled by guard cell exchange)
        {
            Box bx = mf[mfi].box();
            const Box realspace_bx = bx.enclosedCells(); // discards last point in each nodal direction
            AMREX_ALWAYS_ASSERT( realspace_bx == tmpRealField[mfi].box() );
            Array4<Real> mf_arr = mf[mfi].array();
            Array4<const Complex> tmp_arr = tmpRealField[mfi].array();
            ParallelFor( realspace_bx,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                mf_arr(i,j,k) = tmp_arr(i,j,k).real();
            });
        }
    }
}


SpectralField&
SpectralFieldData::getSpectralField( const int field_index )
{
    switch(field_index)
    {
        case SpectralFieldIndex::Ex : return Ex;
        case SpectralFieldIndex::Ey : return Ey;
        case SpectralFieldIndex::Ez : return Ez;
        case SpectralFieldIndex::Bx : return Bx;
        case SpectralFieldIndex::By : return By;
        case SpectralFieldIndex::Bz : return Bz;
        case SpectralFieldIndex::Jx : return Jx;
        case SpectralFieldIndex::Jy : return Jy;
        case SpectralFieldIndex::Jz : return Jz;
        case SpectralFieldIndex::rho_old : return rho_old;
        case SpectralFieldIndex::rho_new : return rho_new;
        default : return tmpSpectralField; // For synthax; should not occur in practice
    }
}
