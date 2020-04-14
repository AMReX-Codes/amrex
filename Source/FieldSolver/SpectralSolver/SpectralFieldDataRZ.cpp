/* Copyright 2019-2020 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "SpectralFieldDataRZ.H"

#include "WarpX.H"

using amrex::operator""_rt;

/* \brief Initialize fields in spectral space, and FFT plans
 *
 * \param realspace_ba Box array that corresponds to the decomposition
 *  * of the fields in real space (cell-centered ; includes guard cells only in z)
 * \param k_space Defined the domain of the k space
 * \param dm Indicates which MPI proc owns which box, in realspace_ba
 * \param n_field_required Specifies the number of fields that will be transformed
 * \param n_modes Number of cylindrical modes
 * */
SpectralFieldDataRZ::SpectralFieldDataRZ (amrex::BoxArray const & realspace_ba,
                                          SpectralKSpaceRZ const & k_space,
                                          amrex::DistributionMapping const & dm,
                                          int const n_field_required,
                                          int const n_modes,
                                          int const lev)
    : n_rz_azimuthal_modes(n_modes)
{
    amrex::BoxArray const & spectralspace_ba = k_space.spectralspace_ba;

    // Allocate the arrays that contain the fields in spectral space.
    // SpectralField is comparable to a MultiFab but stores complex numbers.
    // This stores all of the transformed fields in one place, with the last dimension
    // being the list of fields, defined by SpectralFieldIndex, for all of the modes.
    // The fields of each mode are grouped together, so that the index of a
    // field for a specific mode is given by field_index + mode*n_fields.
    fields = SpectralField(spectralspace_ba, dm, n_rz_azimuthal_modes*n_field_required, 0);

    // Allocate temporary arrays - in real space and spectral space.
    // These complex arrays will store the data just before/after the z FFT.
    // Note that the realspace_ba should not include the radial guard cells.
    tempHTransformed = SpectralField(realspace_ba, dm, n_rz_azimuthal_modes, 0);
    tmpSpectralField = SpectralField(spectralspace_ba, dm, n_rz_azimuthal_modes, 0);

    // By default, we assume the z FFT is done from/to a nodal grid in real space.
    // It the FFT is performed from/to a cell-centered grid in real space,
    // a correcting "shift" factor must be applied in spectral space.
    zshift_FFTfromCell = k_space.getSpectralShiftFactor(dm, 1,
                                    ShiftType::TransformFromCellCentered);
    zshift_FFTtoCell = k_space.getSpectralShiftFactor(dm, 1,
                                    ShiftType::TransformToCellCentered);

    // Allocate and initialize the FFT plans and Hankel transformer.
    forward_plan = FFTplans(spectralspace_ba, dm);
#ifndef AMREX_USE_GPU
    // The backward plan is not needed with GPU since it would be the same
    // as the forward plan anyway.
    backward_plan = FFTplans(spectralspace_ba, dm);
#endif
    multi_spectral_hankel_transformer = MultiSpectralHankelTransformer(spectralspace_ba, dm);

    // Loop over boxes and allocate the corresponding plan
    // for each box owned by the local MPI proc.
    for (amrex::MFIter mfi(spectralspace_ba, dm); mfi.isValid(); ++mfi){
        amrex::IntVect grid_size = realspace_ba[mfi].length();
#ifdef AMREX_USE_GPU
        // Create cuFFT plan.
        // This is alway complex to complex.
        // This plan is for one azimuthal mode only.
        cufftResult result;
        int fft_length[] = {grid_size[1]};
        int inembed[] = {grid_size[1]};
        int istride = grid_size[0];
        int idist = 1;
        int onembed[] = {grid_size[1]};
        int ostride = grid_size[0];
        int odist = 1;
        int batch = grid_size[0]; // number of ffts
        result = cufftPlanMany(&forward_plan[mfi], 1, fft_length, inembed, istride, idist,
                               onembed, ostride, odist, CUFFT_Z2Z, batch);
        if (result != CUFFT_SUCCESS) {
           amrex::Print() << " cufftPlanMany failed! \n";
        }
        // The backward plane is the same as the forward since the direction is passed when executed.

#else
        // Create FFTW plans.
        fftw_iodim dims[1];
        fftw_iodim howmany_dims[2];
        dims[0].n = grid_size[1];
        dims[0].is = grid_size[0];
        dims[0].os = grid_size[0];
        howmany_dims[0].n = n_rz_azimuthal_modes;
        howmany_dims[0].is = grid_size[0]*grid_size[1];
        howmany_dims[0].os = grid_size[0]*grid_size[1];
        howmany_dims[1].n = grid_size[0];
        howmany_dims[1].is = 1;
        howmany_dims[1].os = 1;
        forward_plan[mfi] =
            // Note that AMReX FAB are Fortran-order.
            fftw_plan_guru_dft(1, // int rank
                               dims,
                               2, // int howmany_rank,
                               howmany_dims,
                               reinterpret_cast<fftw_complex*>(tempHTransformed[mfi].dataPtr()), // fftw_complex *in
                               reinterpret_cast<fftw_complex*>(tmpSpectralField[mfi].dataPtr()), // fftw_complex *out
                               FFTW_FORWARD, // int sign
                               FFTW_ESTIMATE); // unsigned flags
        backward_plan[mfi] =
            fftw_plan_guru_dft(1, // int rank
                               dims,
                               2, // int howmany_rank,
                               howmany_dims,
                               reinterpret_cast<fftw_complex*>(tmpSpectralField[mfi].dataPtr()), // fftw_complex *in
                               reinterpret_cast<fftw_complex*>(tempHTransformed[mfi].dataPtr()), // fftw_complex *out
                               FFTW_BACKWARD, // int sign
                               FFTW_ESTIMATE); // unsigned flags
#endif

        // Create the Hankel transformer for each box.
        std::array<amrex::Real,3> xmax = WarpX::UpperCorner(mfi.tilebox(), lev);
        multi_spectral_hankel_transformer[mfi] = SpectralHankelTransformer(grid_size[0], n_rz_azimuthal_modes, xmax[0]);
    }
}


SpectralFieldDataRZ::~SpectralFieldDataRZ()
{
    if (fields.size() > 0){
        for (amrex::MFIter mfi(fields); mfi.isValid(); ++mfi){
#ifdef AMREX_USE_GPU
            // Destroy cuFFT plans.
            cufftDestroy(forward_plan[mfi]);
            // cufftDestroy(backward_plan[mfi]); // This was never allocated.
#else
            // Destroy FFTW plans.
            fftw_destroy_plan(forward_plan[mfi]);
            fftw_destroy_plan(backward_plan[mfi]);
#endif
        }
    }
}

/* \brief Z Transform the FAB to spectral space,
 *  and store the corresponding result internally
 *  (in the spectral field specified by `field_index`)
 *  The input, tempHTransformedSplit, is the complex, Hankel transformed
 *  data, which is stored wih the real and imaginary parts split.
 *  The input should include the imaginary component of mode 0
 *  (even though it is all zeros). */
void
SpectralFieldDataRZ::FABZForwardTransform (amrex::MFIter const & mfi,
                                           amrex::MultiFab const & tempHTransformedSplit,
                                           int const field_index, const bool is_nodal_z)
{
    // Copy the split complex to the interleaved complex.

    amrex::Box const& realspace_bx = tempHTransformed[mfi].box();

    amrex::Array4<const amrex::Real> const& split_arr = tempHTransformedSplit[mfi].array();
    amrex::Array4<Complex> const& complex_arr = tempHTransformed[mfi].array();

    int const modes = n_rz_azimuthal_modes;
    ParallelFor(realspace_bx, modes,
    [=] AMREX_GPU_DEVICE(int i, int j, int k, int mode) noexcept {
        int const mode_r = 2*mode;
        int const mode_i = 2*mode + 1;
        complex_arr(i,j,k,mode) = Complex{split_arr(i,j,k,mode_r), split_arr(i,j,k,mode_i)};
    });

    // Perform Fourier transform from `tempHTransformed` to `tmpSpectralField`.
#ifdef AMREX_USE_GPU
    // Perform Fast Fourier Transform on GPU using cuFFT.
    // Make sure that this is done on the same
    // GPU stream as the above copy.
    cufftResult result;
    cudaStream_t stream = amrex::Gpu::Device::cudaStream();
    cufftSetStream(forward_plan[mfi], stream);
    for (int mode=0 ; mode < n_rz_azimuthal_modes ; mode++) {
        result = cufftExecZ2Z(forward_plan[mfi],
                              reinterpret_cast<cuDoubleComplex*>(tempHTransformed[mfi].dataPtr(mode)), // cuDoubleComplex *in
                              reinterpret_cast<cuDoubleComplex*>(tmpSpectralField[mfi].dataPtr(mode)), // cuDoubleComplex *out
                              CUFFT_FORWARD);
        if (result != CUFFT_SUCCESS) {
           amrex::Print() << " forward transform using cufftExecZ2Z failed ! \n";
        }
    }
#else
    fftw_execute(forward_plan[mfi]);
#endif

    // Copy the spectral-space field `tmpSpectralField` to the appropriate
    // index of the FabArray `fields` (specified by `field_index`)
    // and apply correcting shift factor if the real space data comes
    // from a cell-centered grid in real space instead of a nodal grid.
    amrex::Array4<const Complex> const& tmp_arr = tmpSpectralField[mfi].array();
    amrex::Array4<Complex> const& fields_arr = fields[mfi].array();
    Complex const* zshift_arr = zshift_FFTfromCell[mfi].dataPtr();

    // Loop over indices within one box, all components.
    // The fields are organized so that the fields for each mode
    // are grouped together in memory.
    amrex::Box const& spectralspace_bx = tmpSpectralField[mfi].box();
    int const nz = spectralspace_bx.length(1);
    amrex::Real inv_nz = 1._rt/nz;
    constexpr int n_fields = SpectralFieldIndex::n_fields;

    ParallelFor(spectralspace_bx, modes,
    [=] AMREX_GPU_DEVICE(int i, int j, int k, int mode) noexcept {
        Complex spectral_field_value = tmp_arr(i,j,k,mode);
        // Apply proper shift.
        if (is_nodal_z==false) spectral_field_value *= zshift_arr[j];
        // Copy field into the correct index.
        int const ic = field_index + mode*n_fields;
        fields_arr(i,j,k,ic) = spectral_field_value*inv_nz;
    });
}

/* \brief Backward Z Transform the data from the fields
 * (in the spectral field specified by `field_index`)
 * to physical space, and return the resulting FArrayBox.
 *  The output, tempHTransformedSplit, is the complex, Hankel transformed
 *  data, which is stored wih the real and imaginary parts split.
 *  The output includes the imaginary component of mode 0
 *  (even though it is all zeros). */
void
SpectralFieldDataRZ::FABZBackwardTransform (amrex::MFIter const & mfi, int const field_index,
                                            amrex::MultiFab & tempHTransformedSplit,
                                            const bool is_nodal_z)
{
    // Copy the spectral-space field from the appropriate index of the FabArray
    // `fields` (specified by `field_index`) to field `tmpSpectralField`
    // and apply correcting shift factor if the real space data is on
    // a cell-centered grid in real space instead of a nodal grid.
    amrex::Array4<const Complex> const& fields_arr = fields[mfi].array();
    amrex::Array4<Complex> const& tmp_arr = tmpSpectralField[mfi].array();
    Complex const* zshift_arr = zshift_FFTtoCell[mfi].dataPtr();

    // Loop over indices within one box, all components.
    amrex::Box const& spectralspace_bx = tmpSpectralField[mfi].box();

    int const modes = n_rz_azimuthal_modes;
    constexpr int n_fields = SpectralFieldIndex::n_fields;
    ParallelFor(spectralspace_bx, modes,
    [=] AMREX_GPU_DEVICE(int i, int j, int k, int mode) noexcept {
        int const ic = field_index + mode*n_fields;
        Complex spectral_field_value = fields_arr(i,j,k,ic);
        // Apply proper shift.
        if (is_nodal_z==false) spectral_field_value *= zshift_arr[j];
        // Copy field into the right index.
        tmp_arr(i,j,k,mode) = spectral_field_value;
    });

    // Perform Fourier transform from `tmpSpectralField` to `tempHTransformed`.
#ifdef AMREX_USE_GPU
    // Perform Fast Fourier Transform on GPU using cuFFT.
    // Make sure that this is done on the same
    // GPU stream as the above copy.
    cufftResult result;
    cudaStream_t stream = amrex::Gpu::Device::cudaStream();
    cufftSetStream(forward_plan[mfi], stream);
    for (int mode=0 ; mode < n_rz_azimuthal_modes ; mode++) {
        result = cufftExecZ2Z(forward_plan[mfi],
                              reinterpret_cast<cuDoubleComplex*>(tmpSpectralField[mfi].dataPtr(mode)), // cuDoubleComplex *in
                              reinterpret_cast<cuDoubleComplex*>(tempHTransformed[mfi].dataPtr(mode)), // cuDoubleComplex *out
                              CUFFT_INVERSE);
        if (result != CUFFT_SUCCESS) {
           amrex::Print() << " backwardtransform using cufftExecZ2Z failed ! \n";
        }
    }
#else
    fftw_execute(backward_plan[mfi]);
#endif

    // Copy the interleaved complex to the split complex.
    amrex::Box const& realspace_bx = tempHTransformed[mfi].box();

    amrex::Array4<amrex::Real> const& split_arr = tempHTransformedSplit[mfi].array();
    amrex::Array4<const Complex> const& complex_arr = tempHTransformed[mfi].array();

    ParallelFor(realspace_bx, modes,
    [=] AMREX_GPU_DEVICE(int i, int j, int k, int mode) noexcept {
        int const mode_r = 2*mode;
        int const mode_i = 2*mode + 1;
        split_arr(i,j,k,mode_r) = complex_arr(i,j,k,mode).real();
        split_arr(i,j,k,mode_i) = complex_arr(i,j,k,mode).imag();
    });

}

/* \brief Transform the component `i_comp` of MultiFab `field_mf`
 *  to spectral space, and store the corresponding result internally
 *  (in the spectral field specified by `field_index`) */
void
SpectralFieldDataRZ::ForwardTransform (amrex::MultiFab const & field_mf, int const field_index,
                                       int const i_comp)
{
    // Check field index type, in order to apply proper shift in spectral space.
    // Only cell centered in r is supported.
    bool const is_nodal_z = field_mf.is_nodal(1);

    int const ncomp = 2*n_rz_azimuthal_modes - 1;

    // This will hold the Hankel transformed data, with the real and imaginary parts split.
    // A full multifab is created so that each GPU stream has its own temp space.
    amrex::MultiFab tempHTransformedSplit(tempHTransformed.boxArray(), tempHTransformed.DistributionMap(), 2*n_rz_azimuthal_modes, 0);

    // Loop over boxes.
    for (amrex::MFIter mfi(field_mf); mfi.isValid(); ++mfi){

        // Perform the Hankel transform first.
        // tempHTransformedSplit includes the imaginary component of mode 0.
        // field_mf does not.
        amrex::Box const& realspace_bx = tempHTransformed[mfi].box();
        amrex::FArrayBox field_comp(field_mf[mfi], amrex::make_alias, i_comp*ncomp, ncomp);
        multi_spectral_hankel_transformer[mfi].PhysicalToSpectral_Scalar(realspace_bx, field_comp, tempHTransformedSplit[mfi]);

        FABZForwardTransform(mfi, tempHTransformedSplit, field_index, is_nodal_z);

    }
}

/* \brief Transform the coupled components of MultiFabs `field_mf_r` and `field_mf_t`
 *  to spectral space, and store the corresponding result internally
 *  (in the spectral fields specified by `field_index_r` and `field_index_t`) */
void
SpectralFieldDataRZ::ForwardTransform (amrex::MultiFab const & field_mf_r, int const field_index_r,
                                       amrex::MultiFab const & field_mf_t, int const field_index_t)
{
    // Check field index type, in order to apply proper shift in spectral space.
    // Only cell centered in r is supported.
    bool const is_nodal_z = field_mf_r.is_nodal(1);

    // Create copies of the input multifabs. The copies will include the imaginary part of mode 0.
    // Also note that the Hankel transform will overwrite the copies.
    // Full multifabs are created for the temps so that each GPU stream has its own temp space.
    amrex::MultiFab field_mf_r_copy(field_mf_r.boxArray(), field_mf_r.DistributionMap(), 2*n_rz_azimuthal_modes, field_mf_r.nGrowVect());
    amrex::MultiFab field_mf_t_copy(field_mf_t.boxArray(), field_mf_t.DistributionMap(), 2*n_rz_azimuthal_modes, field_mf_t.nGrowVect());
    amrex::MultiFab::Copy(field_mf_r_copy, field_mf_r, 0, 0, 1, field_mf_r.nGrowVect()); // Real part of mode 0
    amrex::MultiFab::Copy(field_mf_t_copy, field_mf_t, 0, 0, 1, field_mf_t.nGrowVect()); // Real part of mode 0
    field_mf_r_copy.setVal(0._rt, 1, 1, field_mf_r.nGrowVect()); // Imaginary part of mode 0
    field_mf_t_copy.setVal(0._rt, 1, 1, field_mf_t.nGrowVect()); // Imaginary part of mode 0
    amrex::MultiFab::Copy(field_mf_r_copy, field_mf_r, 1, 2, 2*n_rz_azimuthal_modes-2, field_mf_r.nGrowVect());
    amrex::MultiFab::Copy(field_mf_t_copy, field_mf_t, 1, 2, 2*n_rz_azimuthal_modes-2, field_mf_t.nGrowVect());

    amrex::MultiFab tempHTransformedSplit_p(tempHTransformed.boxArray(), tempHTransformed.DistributionMap(), 2*n_rz_azimuthal_modes, 0);
    amrex::MultiFab tempHTransformedSplit_m(tempHTransformed.boxArray(), tempHTransformed.DistributionMap(), 2*n_rz_azimuthal_modes, 0);

    // Loop over boxes.
    for (amrex::MFIter mfi(field_mf_r); mfi.isValid(); ++mfi){

        // Perform the Hankel transform first.
        amrex::Box const& realspace_bx = tempHTransformed[mfi].box();
        multi_spectral_hankel_transformer[mfi].PhysicalToSpectral_Vector(realspace_bx,
                                                           field_mf_r_copy[mfi], field_mf_t_copy[mfi],
                                                           tempHTransformedSplit_p[mfi], tempHTransformedSplit_m[mfi]);

        FABZForwardTransform(mfi, tempHTransformedSplit_p, field_index_r, is_nodal_z);
        FABZForwardTransform(mfi, tempHTransformedSplit_m, field_index_t, is_nodal_z);

    }
}

/* \brief Transform spectral field specified by `field_index` back to
 * real space, and store it in the component `i_comp` of `field_mf` */
void
SpectralFieldDataRZ::BackwardTransform (amrex::MultiFab& field_mf, int const field_index,
                                        int const i_comp)
{
    // Check field index type, in order to apply proper shift in spectral space.
    bool const is_nodal_z = field_mf.is_nodal(1);

    int const ncomp = 2*n_rz_azimuthal_modes - 1;

    // A full multifab is created so that each GPU stream has its own temp space.
    amrex::MultiFab tempHTransformedSplit(tempHTransformed.boxArray(), tempHTransformed.DistributionMap(), 2*n_rz_azimuthal_modes, 0);

    // Loop over boxes.
    for (amrex::MFIter mfi(field_mf); mfi.isValid(); ++mfi){

        FABZBackwardTransform(mfi, field_index, tempHTransformedSplit, is_nodal_z);

        // Perform the Hankel inverse transform last.
        // tempHTransformedSplit includes the imaginary component of mode 0.
        // field_mf does not.
        amrex::Box const& realspace_bx = tempHTransformed[mfi].box();
        amrex::FArrayBox field_comp(field_mf[mfi], amrex::make_alias, i_comp*ncomp, ncomp);
        multi_spectral_hankel_transformer[mfi].SpectralToPhysical_Scalar(realspace_bx, tempHTransformedSplit[mfi], field_comp);

    }
}

/* \brief Transform spectral fields specified by `field_index_r` and
 * `field_index_t` back to real space, and store them in `field_mf_r` and `field_mf_t` */
void
SpectralFieldDataRZ::BackwardTransform (amrex::MultiFab& field_mf_r, int const field_index_r,
                                        amrex::MultiFab& field_mf_t, int const field_index_t)
{
    // Check field index type, in order to apply proper shift in spectral space.
    bool const is_nodal_z = field_mf_r.is_nodal(1);

    // Full multifabs are created for the temps so that each GPU stream has its own temp space.
    amrex::MultiFab tempHTransformedSplit_p(tempHTransformed.boxArray(), tempHTransformed.DistributionMap(), 2*n_rz_azimuthal_modes, 0);
    amrex::MultiFab tempHTransformedSplit_m(tempHTransformed.boxArray(), tempHTransformed.DistributionMap(), 2*n_rz_azimuthal_modes, 0);

    // Create copies of the input multifabs. The copies will include the imaginary part of mode 0.
    amrex::MultiFab field_mf_r_copy(field_mf_r.boxArray(), field_mf_r.DistributionMap(), 2*n_rz_azimuthal_modes, field_mf_r.nGrowVect());
    amrex::MultiFab field_mf_t_copy(field_mf_t.boxArray(), field_mf_t.DistributionMap(), 2*n_rz_azimuthal_modes, field_mf_t.nGrowVect());

    // Loop over boxes.
    for (amrex::MFIter mfi(field_mf_r); mfi.isValid(); ++mfi){

        FABZBackwardTransform(mfi, field_index_r, tempHTransformedSplit_p, is_nodal_z);
        FABZBackwardTransform(mfi, field_index_t, tempHTransformedSplit_m, is_nodal_z);

        // Perform the Hankel inverse transform last.
        // tempHTransformedSplit includes the imaginary component of mode 0.
        // field_mf_[ri] do not.
        amrex::Box const& realspace_bx = tempHTransformed[mfi].box();
        multi_spectral_hankel_transformer[mfi].SpectralToPhysical_Vector(realspace_bx,
                                                           tempHTransformedSplit_p[mfi], tempHTransformedSplit_m[mfi],
                                                           field_mf_r_copy[mfi], field_mf_t_copy[mfi]);

        amrex::Array4<amrex::Real> const & field_mf_r_array = field_mf_r[mfi].array();
        amrex::Array4<amrex::Real> const & field_mf_t_array = field_mf_t[mfi].array();
        amrex::Array4<amrex::Real> const & field_mf_r_copy_array = field_mf_r_copy[mfi].array();
        amrex::Array4<amrex::Real> const & field_mf_t_copy_array = field_mf_t_copy[mfi].array();

        ParallelFor(realspace_bx, 2*n_rz_azimuthal_modes-1,
        [=] AMREX_GPU_DEVICE(int i, int j, int k, int icomp) noexcept {
            if (icomp == 0) {
                // mode 0
                field_mf_r_array(i,j,k,icomp) = field_mf_r_copy_array(i,j,k,icomp);
                field_mf_t_array(i,j,k,icomp) = field_mf_t_copy_array(i,j,k,icomp);
            } else {
                field_mf_r_array(i,j,k,icomp) = field_mf_r_copy_array(i,j,k,icomp+1);
                field_mf_t_array(i,j,k,icomp) = field_mf_t_copy_array(i,j,k,icomp+1);
            }
        });

    }

}
