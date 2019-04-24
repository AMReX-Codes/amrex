#include <SpectralFieldData.H>

using namespace amrex;

/* \brief Initialize fields in spectral space, and FFT plans */
SpectralFieldData::SpectralFieldData( const BoxArray& realspace_ba,
                            const SpectralKSpace& k_space,
                            const DistributionMapping& dm )
{
    const BoxArray& spectralspace_ba = k_space.spectralspace_ba;

    // Allocate the arrays that contain the fields in spectral space
    // (one component per field)
    fields = SpectralField(spectralspace_ba, dm,
                            SpectralFieldIndex::n_fields, 0);

    // Allocate temporary arrays - in real space and spectral space
    // These arrays will store the data just before/after the FFT
    tmpRealField = MultiFab(realspace_ba, dm, 1, 0);
    tmpSpectralField = SpectralField(spectralspace_ba, dm, 1, 0);

    // By default, we assume the FFT is done from/to a nodal grid in real space
    // It the FFT is performed from/to a cell-centered grid in real space,
    // a correcting "shift" factor must be applied in spectral space.
    xshift_FFTfromCell = k_space.getSpectralShiftFactor(dm, 0,
                                    ShiftType::TransformFromCellCentered);
    xshift_FFTtoCell = k_space.getSpectralShiftFactor(dm, 0,
                                    ShiftType::TransformToCellCentered);
#if (AMREX_SPACEDIM == 3)
    yshift_FFTfromCell = k_space.getSpectralShiftFactor(dm, 1,
                                    ShiftType::TransformFromCellCentered);
    yshift_FFTtoCell = k_space.getSpectralShiftFactor(dm, 1,
                                    ShiftType::TransformToCellCentered);
    zshift_FFTfromCell = k_space.getSpectralShiftFactor(dm, 2,
                                    ShiftType::TransformFromCellCentered);
    zshift_FFTtoCell = k_space.getSpectralShiftFactor(dm, 2,
                                    ShiftType::TransformToCellCentered);
#else
    zshift_FFTfromCell = k_space.getSpectralShiftFactor(dm, 1,
                                    ShiftType::TransformFromCellCentered);
    zshift_FFTtoCell = k_space.getSpectralShiftFactor(dm, 1,
                                    ShiftType::TransformToCellCentered);
#endif

    // Allocate and initialize the FFT plans
    forward_plan = FFTplans(spectralspace_ba, dm);
    backward_plan = FFTplans(spectralspace_ba, dm);
    // Loop over boxes and allocate the corresponding plan
    // for each box owned by the local MPI proc
    for ( MFIter mfi(spectralspace_ba, dm); mfi.isValid(); ++mfi ){
        // Note: the size of the real-space box and spectral-space box
        // differ when using real-to-complex FFT. When initializing
        // the FFT plan, the valid dimensions are those of the real-space box.
        IntVect fft_size = realspace_ba[mfi].length();
#ifdef AMREX_USE_GPU
        // Add cuFFT-specific code
#else
        // Create FFTW plans
        forward_plan[mfi] =
            // Swap dimensions: AMReX FAB are Fortran-order but FFTW is C-order
#if (AMREX_SPACEDIM == 3)
            fftw_plan_dft_r2c_3d( fft_size[2], fft_size[1], fft_size[0],
#else
            fftw_plan_dft_r2c_2d( fft_size[1], fft_size[0],
#endif
            tmpRealField[mfi].dataPtr(),
            reinterpret_cast<fftw_complex*>( tmpSpectralField[mfi].dataPtr() ),
            FFTW_ESTIMATE );
        backward_plan[mfi] =
            // Swap dimensions: AMReX FAB are Fortran-order but FFTW is C-order
#if (AMREX_SPACEDIM == 3)
            fftw_plan_dft_c2r_3d( fft_size[2], fft_size[1], fft_size[0],
#else
            fftw_plan_dft_c2r_2d( fft_size[1], fft_size[0],
#endif
            reinterpret_cast<fftw_complex*>( tmpSpectralField[mfi].dataPtr() ),
            tmpRealField[mfi].dataPtr(),
            FFTW_ESTIMATE );
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

/* \brief Transform the component `i_comp` of MultiFab `mf`
 *  to spectral space, and store the corresponding result internally
 *  (in the spectral field specified by `field_index`) */
void
SpectralFieldData::ForwardTransform( const MultiFab& mf,
                                     const int field_index,
                                     const int i_comp )
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
            Box realspace_bx = mf[mfi].box(); // Copy the box
            realspace_bx.enclosedCells(); // Discard last point in nodal direction
            AMREX_ALWAYS_ASSERT( realspace_bx == tmpRealField[mfi].box() );
            Array4<const Real> mf_arr = mf[mfi].array();
            Array4<Real> tmp_arr = tmpRealField[mfi].array();
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

        // Copy the spectral-space field `tmpSpectralField` to the appropriate
        // index of the FabArray `fields` (specified by `field_index`)
        // and apply correcting shift factor if the real space data comes
        // from a cell-centered grid in real space instead of a nodal grid.
        {
            Array4<Complex> fields_arr = SpectralFieldData::fields[mfi].array();
            Array4<const Complex> tmp_arr = tmpSpectralField[mfi].array();
            const Complex* xshift_arr = xshift_FFTfromCell[mfi].dataPtr();
#if (AMREX_SPACEDIM == 3)
            const Complex* yshift_arr = yshift_FFTfromCell[mfi].dataPtr();
#endif
            const Complex* zshift_arr = zshift_FFTfromCell[mfi].dataPtr();
            // Loop over indices within one box
            const Box spectralspace_bx = tmpSpectralField[mfi].box();
            ParallelFor( spectralspace_bx,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                Complex spectral_field_value = tmp_arr(i,j,k);
                // Apply proper shift in each dimension
                if (is_nodal_x==false) spectral_field_value *= xshift_arr[i];
#if (AMREX_SPACEDIM == 3)
                if (is_nodal_y==false) spectral_field_value *= yshift_arr[j];
                if (is_nodal_z==false) spectral_field_value *= zshift_arr[k];
#elif (AMREX_SPACEDIM == 2)
                if (is_nodal_z==false) spectral_field_value *= zshift_arr[j];
#endif
                // Copy field into the right index
                fields_arr(i,j,k,field_index) = spectral_field_value;
            });
        }
    }
}


/* \brief Transform spectral field specified by `field_index` back to
 * real space, and store it in the component `i_comp` of `mf` */
void
SpectralFieldData::BackwardTransform( MultiFab& mf,
                                      const int field_index,
                                      const int i_comp )
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

        // Copy the spectral-space field `tmpSpectralField` to the appropriate
        // field (specified by the input argument field_index)
        // and apply correcting shift factor if the field is to be transformed
        // to a cell-centered grid in real space instead of a nodal grid.
        {
            Array4<const Complex> field_arr = SpectralFieldData::fields[mfi].array();
            Array4<Complex> tmp_arr = tmpSpectralField[mfi].array();
            const Complex* xshift_arr = xshift_FFTtoCell[mfi].dataPtr();
#if (AMREX_SPACEDIM == 3)
            const Complex* yshift_arr = yshift_FFTtoCell[mfi].dataPtr();
#endif
            const Complex* zshift_arr = zshift_FFTtoCell[mfi].dataPtr();
            // Loop over indices within one box
            const Box spectralspace_bx = tmpSpectralField[mfi].box();
            ParallelFor( spectralspace_bx,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                Complex spectral_field_value = field_arr(i,j,k,field_index);
                // Apply proper shift in each dimension
                if (is_nodal_x==false) spectral_field_value *= xshift_arr[i];
#if (AMREX_SPACEDIM == 3)
                if (is_nodal_y==false) spectral_field_value *= yshift_arr[j];
                if (is_nodal_z==false) spectral_field_value *= zshift_arr[k];
#elif (AMREX_SPACEDIM == 2)
                if (is_nodal_z==false) spectral_field_value *= zshift_arr[j];
#endif
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

        // Normalize (divide by 1/N) since the FFT+IFFT results in a factor N
        {
            const Box realspace_bx = tmpRealField[mfi].box();
            Array4<Real> mf_arr = mf[mfi].array();
            Array4<const Real> tmp_arr = tmpRealField[mfi].array();
            // Normalization: divide by the number of points in realspace
            const Real inv_N = 1./realspace_bx.numPts();
            ParallelFor( realspace_bx,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                // Copy and normalize field
                mf_arr(i,j,k,i_comp) = inv_N*tmp_arr(i,j,k);
            });
        }
    }
}
