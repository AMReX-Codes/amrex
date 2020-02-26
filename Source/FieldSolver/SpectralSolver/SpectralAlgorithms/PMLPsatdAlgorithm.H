/* Copyright 2019 Axel Huebl, Remi Lehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PML_PSATD_ALGORITHM_H_
#define WARPX_PML_PSATD_ALGORITHM_H_

#include "SpectralBaseAlgorithm.H"

#if WARPX_USE_PSATD

/* \brief Class that updates the field in spectral space
 * and stores the coefficients of the corresponding update equation.
 */
class PMLPsatdAlgorithm : public SpectralBaseAlgorithm
{
    public:
        PMLPsatdAlgorithm(const SpectralKSpace& spectral_kspace,
                         const amrex::DistributionMapping& dm,
                         const int norder_x, const int norder_y,
                         const int norder_z, const bool nodal,
                         const amrex::Real dt);

        void InitializeSpectralCoefficients(
            const SpectralKSpace& spectral_kspace,
            const amrex::DistributionMapping& dm,
            const amrex::Real dt);

        // Redefine functions from base class
        virtual void pushSpectralFields(SpectralFieldData& f) const override final;
        virtual int getRequiredNumberOfFields() const override final {
            return SpectralPMLIndex::n_fields;
        }

    private:
        SpectralRealCoefficients C_coef, S_ck_coef;

};

#endif // WARPX_USE_PSATD
#endif // WARPX_PML_PSATD_ALGORITHM_H_
