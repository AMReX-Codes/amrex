#ifndef WARPX_GALILEAN_ALGORITHM_H_
#define WARPX_GALILEAN_ALGORITHM_H_

#include "SpectralBaseAlgorithm.H"

#if WARPX_USE_PSATD
/* \brief Class that updates the field in spectral space
 * and stores the coefficients of the corresponding update equation.
 */
class GalileanAlgorithm : public SpectralBaseAlgorithm
{
    public:
        GalileanAlgorithm(const SpectralKSpace& spectral_kspace,
                         const amrex::DistributionMapping& dm,
                         const int norder_x, const int norder_y,
                         const int norder_z, const bool nodal,
                         const amrex::Array<amrex::Real,3>& v_galilean,
                         const amrex::Real dt);
        // Redefine update equation from base class
        virtual void pushSpectralFields(SpectralFieldData& f) const override final;
        virtual int getRequiredNumberOfFields() const override final {
            return SpectralFieldIndex::n_fields;
        };
        void InitializeSpectralCoefficients(const SpectralKSpace& spectral_kspace,
                                    const amrex::DistributionMapping& dm,
                                    const amrex::Array<amrex::Real, 3>& v_galilean,
                                    const amrex::Real dt);

    private:
        SpectralRealCoefficients C_coef, S_ck_coef;
        SpectralComplexCoefficients Theta2_coef, X1_coef, X2_coef, X3_coef, X4_coef;
};
#endif // WARPX_USE_PSATD
#endif // WARPX_GALILEAN_ALGORITHM_H_
