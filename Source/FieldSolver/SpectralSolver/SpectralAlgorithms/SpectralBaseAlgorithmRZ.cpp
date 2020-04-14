/* Copyright 2019 Edoardo Zoni
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "SpectralBaseAlgorithmRZ.H"

#include <cmath>

using namespace amrex;

/**
 * \brief Compute spectral divergence of E
 */
void
SpectralBaseAlgorithmRZ::ComputeSpectralDivE (
    SpectralFieldDataRZ& field_data,
    const std::array<std::unique_ptr<amrex::MultiFab>,3>& Efield,
    amrex::MultiFab& divE )
{
    using amrex::operator""_rt;
    using Idx = SpectralFieldIndex;

    // Forward Fourier transform of E
    field_data.ForwardTransform( *Efield[0], Idx::Ex,
                                 *Efield[1], Idx::Ey );
    field_data.ForwardTransform( *Efield[2], Idx::Ez, 0 );

    // Loop over boxes
    for (MFIter mfi(field_data.fields); mfi.isValid(); ++mfi){

        Box const & bx = field_data.fields[mfi].box();

        // Extract arrays for the fields to be updated
        Array4<Complex> fields = field_data.fields[mfi].array();

        // Extract pointers for the k vectors
        auto const & kr = field_data.getKrArray(mfi);
        Real const * kr_arr = kr.dataPtr();
        Real const * modified_kz_arr = modified_kz_vec[mfi].dataPtr();

        int const nr = bx.length(0);
        int const modes = field_data.n_rz_azimuthal_modes;
        constexpr int n_fields = SpectralFieldIndex::n_fields;

        // Loop over indices within one box
        ParallelFor(bx, modes,
        [=] AMREX_GPU_DEVICE(int i, int j, int k, int mode) noexcept
        {
            int const ic1 = Idx::Ex + mode*n_fields;
            int const ic2 = Idx::Ey + mode*n_fields;
            int const ic3 = Idx::Ez + mode*n_fields;

            // Shortcuts for the components of E
            Complex const Ep = fields(i,j,0,ic1);
            Complex const Em = fields(i,j,0,ic2);
            Complex const Ez = fields(i,j,0,ic3);

            // k vector values
            int const ir = i + nr*mode;
            Real const kr = kr_arr[ir];
            Real const kz = modified_kz_arr[j];
            Complex const I = Complex{0._rt,1._rt};

            // div(E) in Fourier space
            int const ic = Idx::divE + mode*n_fields;
            fields(i,j,0,ic) = kr*(Ep - Em) + I*kz*Ez;
        });
    }

    // Backward Fourier transform
    field_data.BackwardTransform( divE, Idx::divE, 0 );
};
