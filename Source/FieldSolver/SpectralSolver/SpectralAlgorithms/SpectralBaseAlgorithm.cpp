/* Copyright 2019 Edoardo Zoni
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "SpectralBaseAlgorithm.H"
#include <cmath>

using namespace amrex;

/**
 * \brief Compute spectral divergence of E
 */
void
SpectralBaseAlgorithm::ComputeSpectralDivE (
    SpectralFieldData& field_data,
    const std::array<std::unique_ptr<amrex::MultiFab>,3>& Efield,
    amrex::MultiFab& divE )
{
    using Idx = SpectralFieldIndex;

    // Forward Fourier transform of E
    field_data.ForwardTransform( *Efield[0], Idx::Ex, 0 );
    field_data.ForwardTransform( *Efield[1], Idx::Ey, 0 );
    field_data.ForwardTransform( *Efield[2], Idx::Ez, 0 );

    // Loop over boxes
    for (MFIter mfi(field_data.fields); mfi.isValid(); ++mfi){

        const Box& bx = field_data.fields[mfi].box();

        // Extract arrays for the fields to be updated
        Array4<Complex> fields = field_data.fields[mfi].array();
        // Extract pointers for the k vectors
        const Real* modified_kx_arr = modified_kx_vec[mfi].dataPtr();
#if (AMREX_SPACEDIM==3)
        const Real* modified_ky_arr = modified_ky_vec[mfi].dataPtr();
#endif
        const Real* modified_kz_arr = modified_kz_vec[mfi].dataPtr();

        // Loop over indices within one box
        ParallelFor(bx,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            using Idx = SpectralFieldIndex;
            // Shortcuts for the components of E
            const Complex Ex = fields(i,j,k,Idx::Ex);
            const Complex Ey = fields(i,j,k,Idx::Ey);
            const Complex Ez = fields(i,j,k,Idx::Ez);
            // k vector values
            const Real kx = modified_kx_arr[i];
#if (AMREX_SPACEDIM==3)
            const Real ky = modified_ky_arr[j];
            const Real kz = modified_kz_arr[k];
#else
            constexpr Real ky = 0;
            const Real kz = modified_kz_arr[j];
#endif
            const Complex I = Complex{0,1};

            // div(E) in Fourier space
            fields(i,j,k,Idx::divE) = I*(kx*Ex+ky*Ey+kz*Ez);
        });
    }

    // Backward Fourier transform
    field_data.BackwardTransform( divE, Idx::divE, 0 );
};
