#ifndef WARPX_PSATD_ALGORITHM_H_
#define WARPX_PSATD_ALGORITHM_H_

//#include <SpectralKSpace.H>
//#include <SpectralFieldData.H>
#include <SpectralBaseAlgorithm.H>

/* \brief Class that updates the field in spectral space
 * and stores the coefficients of the corresponding update equation.
 */
class PsatdAlgorithm : public SpectralBaseAlgorithm
{
    using SpectralCoefficients = amrex::FabArray< amrex::BaseFab <amrex::Real> >;

    public:
        PsatdAlgorithm(const SpectralKSpace& spectral_kspace,
                         const amrex::DistributionMapping& dm,
                         const int norder_x, const int norder_y,
                         const int norder_z, const bool nodal, const amrex::Real dt);

        void InitializeSpectralCoefficients(const SpectralKSpace& spectral_kspace,
                                    const amrex::DistributionMapping& dm, 
                                    const amrex::Real dt);

        void pushSpectralFields(SpectralFieldData& f) const override final;

    private:
        // Modified finite-order vectors
        SpectralCoefficients C_coef, S_ck_coef, X1_coef, X2_coef, X3_coef;
};

#endif // WARPX_PSATD_ALGORITHM_H_
