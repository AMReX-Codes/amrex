/* Copyright 2019-2020 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "PsatdAlgorithmRZ.H"
#include "Utils/WarpXConst.H"

#include <cmath>

using amrex::operator""_rt;


/* \brief Initialize coefficients for the update equation */
PsatdAlgorithmRZ::PsatdAlgorithmRZ (SpectralKSpaceRZ const & spectral_kspace,
                                    amrex::DistributionMapping const & dm,
                                    int const n_rz_azimuthal_modes, int const norder_z,
                                    bool const nodal, amrex::Real const dt_step)
     // Initialize members of base class
     : SpectralBaseAlgorithmRZ(spectral_kspace, dm,
                               norder_z, nodal),
       dt(dt_step)
{

    // Allocate the arrays of coefficients
    amrex::BoxArray const & ba = spectral_kspace.spectralspace_ba;
    C_coef = SpectralCoefficients(ba, dm, n_rz_azimuthal_modes, 0);
    S_ck_coef = SpectralCoefficients(ba, dm, n_rz_azimuthal_modes, 0);
    X1_coef = SpectralCoefficients(ba, dm, n_rz_azimuthal_modes, 0);
    X2_coef = SpectralCoefficients(ba, dm, n_rz_azimuthal_modes, 0);
    X3_coef = SpectralCoefficients(ba, dm, n_rz_azimuthal_modes, 0);

    coefficients_initialized = false;
}

/* Advance the E and B field in spectral space (stored in `f`)
 * over one time step */
void
PsatdAlgorithmRZ::pushSpectralFields(SpectralFieldDataRZ & f)
{

    if (not coefficients_initialized) {
        // This is called from here since it needs the kr values
        // which can be obtained from the SpectralFieldDataRZ
        InitializeSpectralCoefficients(f);
        coefficients_initialized = true;
    }

    // Loop over boxes
    for (amrex::MFIter mfi(f.fields); mfi.isValid(); ++mfi){

        amrex::Box const & bx = f.fields[mfi].box();

        // Extract arrays for the fields to be updated
        amrex::Array4<Complex> const& fields = f.fields[mfi].array();
        // Extract arrays for the coefficients
        amrex::Array4<const amrex::Real> const& C_arr = C_coef[mfi].array();
        amrex::Array4<const amrex::Real> const& S_ck_arr = S_ck_coef[mfi].array();
        amrex::Array4<const amrex::Real> const& X1_arr = X1_coef[mfi].array();
        amrex::Array4<const amrex::Real> const& X2_arr = X2_coef[mfi].array();
        amrex::Array4<const amrex::Real> const& X3_arr = X3_coef[mfi].array();

        // Extract pointers for the k vectors
        auto const & kr_modes = f.getKrArray(mfi);
        amrex::Real const* kr_arr = kr_modes.dataPtr();
        amrex::Real const* modified_kz_arr = modified_kz_vec[mfi].dataPtr();
        int const nr = bx.length(0);

        // Loop over indices within one box
        // Note that k = 0
        int const modes = f.n_rz_azimuthal_modes;
        amrex::ParallelFor(bx, modes,
        [=] AMREX_GPU_DEVICE(int i, int j, int k, int mode) noexcept
        {

            // All of the fields of each mode are grouped together
            using Idx = SpectralFieldIndex;
            auto const Ep_m = Idx::Ex + Idx::n_fields*mode;
            auto const Em_m = Idx::Ey + Idx::n_fields*mode;
            auto const Ez_m = Idx::Ez + Idx::n_fields*mode;
            auto const Bp_m = Idx::Bx + Idx::n_fields*mode;
            auto const Bm_m = Idx::By + Idx::n_fields*mode;
            auto const Bz_m = Idx::Bz + Idx::n_fields*mode;
            auto const Jp_m = Idx::Jx + Idx::n_fields*mode;
            auto const Jm_m = Idx::Jy + Idx::n_fields*mode;
            auto const Jz_m = Idx::Jz + Idx::n_fields*mode;
            auto const rho_old_m = Idx::rho_old + Idx::n_fields*mode;
            auto const rho_new_m = Idx::rho_new + Idx::n_fields*mode;

            // Record old values of the fields to be updated
            Complex const Ep_old = fields(i,j,k,Ep_m);
            Complex const Em_old = fields(i,j,k,Em_m);
            Complex const Ez_old = fields(i,j,k,Ez_m);
            Complex const Bp_old = fields(i,j,k,Bp_m);
            Complex const Bm_old = fields(i,j,k,Bm_m);
            Complex const Bz_old = fields(i,j,k,Bz_m);
            // Shortcut for the values of J and rho
            Complex const Jp = fields(i,j,k,Jp_m);
            Complex const Jm = fields(i,j,k,Jm_m);
            Complex const Jz = fields(i,j,k,Jz_m);
            Complex const rho_old = fields(i,j,k,rho_old_m);
            Complex const rho_new = fields(i,j,k,rho_new_m);

            // k vector values, and coefficients
            // The k values for each mode are grouped together
            int const ir = i + nr*mode;
            amrex::Real const kr = kr_arr[ir];
            amrex::Real const kz = modified_kz_arr[j];

            constexpr amrex::Real c2 = PhysConst::c*PhysConst::c;
            constexpr amrex::Real inv_ep0 = 1._rt/PhysConst::ep0;
            Complex const I = Complex{0._rt,1._rt};
            amrex::Real const C = C_arr(i,j,k,mode);
            amrex::Real const S_ck = S_ck_arr(i,j,k,mode);
            amrex::Real const X1 = X1_arr(i,j,k,mode);
            amrex::Real const X2 = X2_arr(i,j,k,mode);
            amrex::Real const X3 = X3_arr(i,j,k,mode);

            // Update E (see WarpX online documentation: theory section)
            fields(i,j,k,Ep_m) = C*Ep_old
                        + S_ck*(-c2*I*kr/2._rt*Bz_old + c2*kz*Bp_old - inv_ep0*Jp)
                        + kr*(X2*rho_new - X3*rho_old);
            fields(i,j,k,Em_m) = C*Em_old
                        + S_ck*(-c2*I*kr/2._rt*Bz_old - c2*kz*Bm_old - inv_ep0*Jm)
                        - kr*(X2*rho_new - X3*rho_old);
            fields(i,j,k,Ez_m) = C*Ez_old
                        + S_ck*(c2*I*kr*Bp_old + c2*I*kr*Bm_old - inv_ep0*Jz)
                        - I*kz*(X2*rho_new - X3*rho_old);
            // Update B (see WarpX online documentation: theory section)
            fields(i,j,k,Bp_m) = C*Bp_old
                        - S_ck*(-I*kr/2._rt*Ez_old + kz*Ep_old)
                        + X1*(-I*kr/2._rt*Jz + kz*Jp);
            fields(i,j,k,Bm_m) = C*Bm_old
                        - S_ck*(-I*kr/2._rt*Ez_old - kz*Em_old)
                        + X1*(-I*kr/2._rt*Jz - kz*Jm);
            fields(i,j,k,Bz_m) = C*Bz_old
                        - S_ck*I*(kr*Ep_old + kr*Em_old)
                        + X1*I*(kr*Jp + kr*Jm);
        });
    }
};

void PsatdAlgorithmRZ::InitializeSpectralCoefficients (SpectralFieldDataRZ const & f)
{

    // Fill them with the right values:
    // Loop over boxes and allocate the corresponding coefficients
    // for each box owned by the local MPI proc
    for (amrex::MFIter mfi(f.fields); mfi.isValid(); ++mfi){

        amrex::Box const & bx = f.fields[mfi].box();

        // Extract pointers for the k vectors
        amrex::Real const* const modified_kz = modified_kz_vec[mfi].dataPtr();

        // Extract arrays for the coefficients
        amrex::Array4<amrex::Real> const& C = C_coef[mfi].array();
        amrex::Array4<amrex::Real> const& S_ck = S_ck_coef[mfi].array();
        amrex::Array4<amrex::Real> const& X1 = X1_coef[mfi].array();
        amrex::Array4<amrex::Real> const& X2 = X2_coef[mfi].array();
        amrex::Array4<amrex::Real> const& X3 = X3_coef[mfi].array();

        auto const & kr_modes = f.getKrArray(mfi);
        amrex::Real const* kr_arr = kr_modes.dataPtr();
        int const nr = bx.length(0);
        amrex::Real const dt_temp = dt;

        // Loop over indices within one box
        int const modes = f.n_rz_azimuthal_modes;
        amrex::ParallelFor(bx, modes,
        [=] AMREX_GPU_DEVICE(int i, int j, int k, int mode) noexcept
        {
            // Calculate norm of vector
            int const ir = i + nr*mode;
            amrex::Real const kr = kr_arr[ir];
            amrex::Real const kz = modified_kz[j];
            amrex::Real const k_norm = std::sqrt(kr*kr + kz*kz);

            // Calculate coefficients
            constexpr amrex::Real c = PhysConst::c;
            constexpr amrex::Real ep0 = PhysConst::ep0;
            if (k_norm != 0){
                C(i,j,k,mode) = std::cos(c*k_norm*dt_temp);
                S_ck(i,j,k,mode) = std::sin(c*k_norm*dt_temp)/(c*k_norm);
                X1(i,j,k,mode) = (1._rt - C(i,j,k,mode))/(ep0 * c*c * k_norm*k_norm);
                X2(i,j,k,mode) = (1._rt - S_ck(i,j,k,mode)/dt_temp)/(ep0 * k_norm*k_norm);
                X3(i,j,k,mode) = (C(i,j,k,mode) - S_ck(i,j,k,mode)/dt_temp)/(ep0 * k_norm*k_norm);
            } else { // Handle k_norm = 0, by using the analytical limit
                C(i,j,k,mode) = 1._rt;
                S_ck(i,j,k,mode) = dt_temp;
                X1(i,j,k,mode) = 0.5_rt * dt_temp*dt_temp / ep0;
                X2(i,j,k,mode) = c*c * dt_temp*dt_temp / (6._rt*ep0);
                X3(i,j,k,mode) = - c*c * dt_temp*dt_temp / (3._rt*ep0);
            }
        });
     }
}
