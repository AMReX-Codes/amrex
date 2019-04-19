#ifndef WARPX_PSATD_ALGORITHM_H_
#define WARPX_PSATD_ALGORITHM_H_

#include <SpectralKSpace.H>
#include <SpectralFieldData.H>

/* TODO: Write documentation
*/
class PsatdAlgorithm
{
    using SpectralCoefficients = amrex::FabArray< amrex::BaseFab <amrex::Real> >;

    public:
        PsatdAlgorithm( const SpectralKSpace& spectral_kspace,
                         const amrex::DistributionMapping& dm,
                         const int norder_x, const int norder_y,
                         const int norder_z, const amrex::Real dt);
        PsatdAlgorithm() = default; // Default constructor
        PsatdAlgorithm& operator=(PsatdAlgorithm&& algorithm) = default; // Default move assignment
        void pushSpectralFields( SpectralFieldData& f ) const;

    private:
        // Modified finite-order vectors
        SpectralKVector modified_kx_vec, modified_ky_vec, modified_kz_vec;
        SpectralCoefficients C_coef, S_ck_coef, X1_coef, X2_coef, X3_coef;
};

#endif // WARPX_PSATD_ALGORITHM_H_
