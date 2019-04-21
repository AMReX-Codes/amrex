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

    // Allocate the coefficients that allow to shift between nodal and cell-centered
    xshift_C2N = k_space.getSpectralShiftFactor(dm, 0, ShiftType::CenteredToNodal);
    xshift_N2C = k_space.getSpectralShiftFactor(dm, 0, ShiftType::NodalToCentered);
#if (AMREX_SPACEDIM == 3)
    yshift_C2N = k_space.getSpectralShiftFactor(dm, 1, ShiftType::CenteredToNodal);
    yshift_N2C = k_space.getSpectralShiftFactor(dm, 1, ShiftType::NodalToCentered);
    zshift_C2N = k_space.getSpectralShiftFactor(dm, 2, ShiftType::CenteredToNodal);
    zshift_N2C = k_space.getSpectralShiftFactor(dm, 2, ShiftType::NodalToCentered);
#else
    zshift_C2N = k_space.getSpectralShiftFactor(dm, 1, ShiftType::CenteredToNodal);
    zshift_N2C = k_space.getSpectralShiftFactor(dm, 1, ShiftType::NodalToCentered);
#endif

    // Allocate and initialize the FFT plans
    forward_plan = FFTplans(spectralspace_ba, dm);
    backward_plan = FFTplans(spectralspace_ba, dm);
    for ( MFIter mfi(spectralspace_ba, dm); mfi.isValid(); ++mfi ){
        Box bx = spectralspace_ba[mfi];
#ifdef AMREX_USE_GPU
        // Add cuFFT-specific code
#else
        // Create FFTW plans
        forward_plan[mfi] =
#if (AMREX_SPACEDIM == 3) // Swap dimensions: AMReX data is Fortran-order, but FFTW is C-order
            fftw_plan_dft_3d( bx.length(2), bx.length(1), bx.length(0),
#else
            fftw_plan_dft_2d( bx.length(1), bx.length(0),
#endif
            reinterpret_cast<fftw_complex*>( tmpRealField[mfi].dataPtr() ),
            reinterpret_cast<fftw_complex*>( tmpSpectralField[mfi].dataPtr() ),
            FFTW_FORWARD, FFTW_ESTIMATE );
        backward_plan[mfi] =
#if (AMREX_SPACEDIM == 3) // Swap dimensions: AMReX data is Fortran-order, but FFTW is C-order
            fftw_plan_dft_3d( bx.length(2), bx.length(1), bx.length(0),
#else
            fftw_plan_dft_2d( bx.length(1), bx.length(0),
#endif
            reinterpret_cast<fftw_complex*>( tmpSpectralField[mfi].dataPtr() ),
            reinterpret_cast<fftw_complex*>( tmpRealField[mfi].dataPtr() ),
            FFTW_BACKWARD, FFTW_ESTIMATE );
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
SpectralFieldData::ForwardTransform( const MultiFab& mf,
                                     const int field_index, const int i_comp )
{
    // Check field index type, in order to apply proper shift in spectral space
    const bool is_nodal_x = mf.is_nodal(0);
#if (AMREX_SPACEDIM == 3)
    const bool is_nodal_y = mf.is_nodal(1);
    const bool is_nodal_z = mf.is_nodal(2);
#else
    const bool is_nodal_z = mf.is_nodal(1);
#endif

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
                tmp_arr(i,j,k) = mf_arr(i,j,k,i_comp);
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
            const Complex* xshift_C2N_arr = xshift_C2N[mfi].dataPtr();
#if (AMREX_SPACEDIM == 3)
            const Complex* yshift_C2N_arr = yshift_C2N[mfi].dataPtr();
#endif
            const Complex* zshift_C2N_arr = zshift_C2N[mfi].dataPtr();
            // Loop over indices within one box
            const Box spectralspace_bx = tmpSpectralField[mfi].box();
            ParallelFor( spectralspace_bx,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                Complex spectral_field_value = tmp_arr(i,j,k);
                // Apply proper shift in each dimension
                if (is_nodal_x==false) spectral_field_value *= xshift_C2N_arr[i];
#if (AMREX_SPACEDIM == 3)
                if (is_nodal_y==false) spectral_field_value *= yshift_C2N_arr[j];
#endif
                if (is_nodal_z==false) spectral_field_value *= zshift_C2N_arr[k];
                // Copy field into temporary array
                field_arr(i,j,k) = spectral_field_value;
            });
        }
    }
}


/* TODO: Documentation
 */
void
SpectralFieldData::BackwardTransform( MultiFab& mf,
                                      const int field_index, const int i_comp )
{
    // Check field index type, in order to apply proper shift in spectral space
    const bool is_nodal_x = mf.is_nodal(0);
#if (AMREX_SPACEDIM == 3)
    const bool is_nodal_y = mf.is_nodal(1);
    const bool is_nodal_z = mf.is_nodal(2);
#else
    const bool is_nodal_z = mf.is_nodal(1);
#endif

    // Loop over boxes
    for ( MFIter mfi(mf); mfi.isValid(); ++mfi ){

        // Copy the appropriate field (specified by the input argument field_index)
        // to the spectral-space field `tmpSpectralField`
        {
            SpectralField& field = getSpectralField( field_index );
            Array4<const Complex> field_arr = field[mfi].array();
            Array4<Complex> tmp_arr = tmpSpectralField[mfi].array();
            const Complex* xshift_N2C_arr = xshift_N2C[mfi].dataPtr();
#if (AMREX_SPACEDIM == 3)
            const Complex* yshift_N2C_arr = yshift_N2C[mfi].dataPtr();
#endif
            const Complex* zshift_N2C_arr = zshift_N2C[mfi].dataPtr();
            // Loop over indices within one box
            const Box spectralspace_bx = tmpSpectralField[mfi].box();
            ParallelFor( spectralspace_bx,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                Complex spectral_field_value = field_arr(i,j,k);
                // Apply proper shift in each dimension
                if (is_nodal_x==false) spectral_field_value *= xshift_N2C_arr[i];
#if (AMREX_SPACEDIM == 3)
                if (is_nodal_y==false) spectral_field_value *= yshift_N2C_arr[j];
#endif
                if (is_nodal_z==false) spectral_field_value *= zshift_N2C_arr[k];
                // Copy field into temporary array
                tmp_arr(i,j,k) = spectral_field_value;
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
        // Normalize (divide by 1/N) since the FFT result in a factor N
        {
            Box bx = mf[mfi].box();
            const Box realspace_bx = bx.enclosedCells();
            // `enclosedells` discards last point in each nodal direction
            AMREX_ALWAYS_ASSERT( realspace_bx == tmpRealField[mfi].box() );
            Array4<Real> mf_arr = mf[mfi].array();
            Array4<const Complex> tmp_arr = tmpRealField[mfi].array();
            // For normalization: divide by the number of points in the box
            const Real inv_N = 1./bx.numPts();
            ParallelFor( realspace_bx,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                mf_arr(i,j,k,i_comp) = inv_N*tmp_arr(i,j,k).real();
            });
        }
    }
}


SpectralField&
SpectralFieldData::getSpectralField( const int field_index )
{
    switch(field_index)
    {
        case SpectralFieldIndex::Ex : return Ex; break;
        case SpectralFieldIndex::Ey : return Ey; break;
        case SpectralFieldIndex::Ez : return Ez; break;
        case SpectralFieldIndex::Bx : return Bx; break;
        case SpectralFieldIndex::By : return By; break;
        case SpectralFieldIndex::Bz : return Bz; break;
        case SpectralFieldIndex::Jx : return Jx; break;
        case SpectralFieldIndex::Jy : return Jy; break;
        case SpectralFieldIndex::Jz : return Jz; break;
        case SpectralFieldIndex::rho_old : return rho_old; break;
        case SpectralFieldIndex::rho_new : return rho_new; break;
        default : return tmpSpectralField; // For synthax; should not occur in practice
    }
}
