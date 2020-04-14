/* Copyright 2019 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PSATD_ALGORITHM_RZ_H_
#define WARPX_PSATD_ALGORITHM_RZ_H_

#include "SpectralBaseAlgorithmRZ.H"

/* \brief Class that updates the field in spectral space
 * and stores the coefficients of the corresponding update equation.
 */
class PsatdAlgorithmRZ : public SpectralBaseAlgorithmRZ
{

    public:
        PsatdAlgorithmRZ(SpectralKSpaceRZ const & spectral_kspace,
                         amrex::DistributionMapping const & dm,
                         int const n_rz_azimuthal_modes, int const norder_z,
                         bool const nodal, amrex::Real const dt_step);
        // Redefine functions from base class
        virtual void pushSpectralFields(SpectralFieldDataRZ & f) override final;
        virtual int getRequiredNumberOfFields() const override final {
            return SpectralFieldIndex::n_fields;
        }

        void InitializeSpectralCoefficients(SpectralFieldDataRZ const & f);

    private:
        bool coefficients_initialized;
        // Note that dt is saved to use in InitializeSpectralCoefficients
        amrex::Real dt;
        SpectralCoefficients C_coef, S_ck_coef, X1_coef, X2_coef, X3_coef;
};

#endif // WARPX_PSATD_ALGORITHM_RZ_H_
