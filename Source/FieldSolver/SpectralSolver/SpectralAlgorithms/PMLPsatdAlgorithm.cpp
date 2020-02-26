/* Copyright 2019 Remi Lehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "PMLPsatdAlgorithm.H"
#include "Utils/WarpXConst.H"

#include <cmath>


#if WARPX_USE_PSATD

using namespace amrex;

/* \brief Initialize coefficients for the update equation */
PMLPsatdAlgorithm::PMLPsatdAlgorithm(
                         const SpectralKSpace& spectral_kspace,
                         const DistributionMapping& dm,
                         const int norder_x, const int norder_y,
                         const int norder_z, const bool nodal, const Real dt)
     // Initialize members of base class
     : SpectralBaseAlgorithm( spectral_kspace, dm,
                              norder_x, norder_y, norder_z, nodal )
{
    const BoxArray& ba = spectral_kspace.spectralspace_ba;

    // Allocate the arrays of coefficients
    C_coef = SpectralRealCoefficients(ba, dm, 1, 0);
    S_ck_coef = SpectralRealCoefficients(ba, dm, 1, 0);

    InitializeSpectralCoefficients(spectral_kspace, dm, dt);
}

/* Advance the E and B field in spectral space (stored in `f`)
 * over one time step */
void
PMLPsatdAlgorithm::pushSpectralFields(SpectralFieldData& f) const{

    // Loop over boxes
    for (MFIter mfi(f.fields); mfi.isValid(); ++mfi){

        const Box& bx = f.fields[mfi].box();

        // Extract arrays for the fields to be updated
        Array4<Complex> fields = f.fields[mfi].array();
        // Extract arrays for the coefficients
        Array4<const Real> C_arr = C_coef[mfi].array();
        Array4<const Real> S_ck_arr = S_ck_coef[mfi].array();
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
            // Record old values of the fields to be updated
            using Idx = SpectralPMLIndex;
            const Complex Ex_old = fields(i,j,k,Idx::Exy) \
                                 + fields(i,j,k,Idx::Exz);
            const Complex Ey_old = fields(i,j,k,Idx::Eyx) \
                                 + fields(i,j,k,Idx::Eyz);
            const Complex Ez_old = fields(i,j,k,Idx::Ezx) \
                                 + fields(i,j,k,Idx::Ezy);
            const Complex Bx_old = fields(i,j,k,Idx::Bxy) \
                                 + fields(i,j,k,Idx::Bxz);
            const Complex By_old = fields(i,j,k,Idx::Byx) \
                                 + fields(i,j,k,Idx::Byz);
            const Complex Bz_old = fields(i,j,k,Idx::Bzx) \
                                 + fields(i,j,k,Idx::Bzy);
            // k vector values, and coefficients
            const Real kx = modified_kx_arr[i];
#if (AMREX_SPACEDIM==3)
            const Real ky = modified_ky_arr[j];
            const Real kz = modified_kz_arr[k];
#else
            constexpr Real ky = 0;
            const Real kz = modified_kz_arr[j];
#endif
            constexpr Real c2 = PhysConst::c*PhysConst::c;
            const Complex I = Complex{0,1};
            const Real C = C_arr(i,j,k);
            const Real S_ck = S_ck_arr(i,j,k);

            // Update E
            fields(i,j,k,Idx::Exy) = C*fields(i,j,k,Idx::Exy) + S_ck*c2*I*ky*Bz_old;
            fields(i,j,k,Idx::Exz) = C*fields(i,j,k,Idx::Exz) - S_ck*c2*I*kz*By_old;
            fields(i,j,k,Idx::Eyz) = C*fields(i,j,k,Idx::Eyz) + S_ck*c2*I*kz*Bx_old;
            fields(i,j,k,Idx::Eyx) = C*fields(i,j,k,Idx::Eyx) - S_ck*c2*I*kx*Bz_old;
            fields(i,j,k,Idx::Ezx) = C*fields(i,j,k,Idx::Ezx) + S_ck*c2*I*kx*By_old;
            fields(i,j,k,Idx::Ezy) = C*fields(i,j,k,Idx::Ezy) - S_ck*c2*I*ky*Bx_old;
            // Update B
            fields(i,j,k,Idx::Bxy) = C*fields(i,j,k,Idx::Bxy) - S_ck*I*ky*Ez_old;
            fields(i,j,k,Idx::Bxz) = C*fields(i,j,k,Idx::Bxz) + S_ck*I*kz*Ey_old;
            fields(i,j,k,Idx::Byz) = C*fields(i,j,k,Idx::Byz) - S_ck*I*kz*Ex_old;
            fields(i,j,k,Idx::Byx) = C*fields(i,j,k,Idx::Byx) + S_ck*I*kx*Ez_old;
            fields(i,j,k,Idx::Bzx) = C*fields(i,j,k,Idx::Bzx) - S_ck*I*kx*Ey_old;
            fields(i,j,k,Idx::Bzy) = C*fields(i,j,k,Idx::Bzy) + S_ck*I*ky*Ex_old;
        });
    }
};

void PMLPsatdAlgorithm::InitializeSpectralCoefficients (
    const SpectralKSpace& spectral_kspace,
    const amrex::DistributionMapping& dm,
    const amrex::Real dt)
{
    const BoxArray& ba = spectral_kspace.spectralspace_ba;
    // Fill them with the right values:
    // Loop over boxes and allocate the corresponding coefficients
    // for each box owned by the local MPI proc
    for (MFIter mfi(ba, dm); mfi.isValid(); ++mfi){

        const Box& bx = ba[mfi];

        // Extract pointers for the k vectors
        const Real* modified_kx = modified_kx_vec[mfi].dataPtr();
#if (AMREX_SPACEDIM==3)
        const Real* modified_ky = modified_ky_vec[mfi].dataPtr();
#endif
        const Real* modified_kz = modified_kz_vec[mfi].dataPtr();
        // Extract arrays for the coefficients
        Array4<Real> C = C_coef[mfi].array();
        Array4<Real> S_ck = S_ck_coef[mfi].array();

        // Loop over indices within one box
        ParallelFor(bx,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            // Calculate norm of vector
            const Real k_norm = std::sqrt(
                std::pow(modified_kx[i], 2) +
#if (AMREX_SPACEDIM==3)
                std::pow(modified_ky[j], 2) +
                std::pow(modified_kz[k], 2));
#else
                std::pow(modified_kz[j], 2));
#endif

            // Calculate coefficients
            constexpr Real c = PhysConst::c;
            if (k_norm != 0){
                C(i,j,k) = std::cos(c*k_norm*dt);
                S_ck(i,j,k) = std::sin(c*k_norm*dt)/(c*k_norm);
            } else { // Handle k_norm = 0, by using the analytical limit
                C(i,j,k) = 1.;
                S_ck(i,j,k) = dt;
            }
        });
    }
};
#endif // WARPX_USE_PSATD
