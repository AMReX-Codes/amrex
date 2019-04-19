#ifndef WARPX_PSATD_ALGORITHM_H_
#define WARPX_PSATD_ALGORITHM_H_

#include <SpectralKSpace.H>
#include <SpectralFieldData.H>

using namespace amrex;
using namespace Gpu;

/* TODO: Write documentation
*/
class PsatdAlgorithm
{
    using SpectralCoefficients = FabArray<BaseFab<Real>>;

    public:
        PsatdAlgorithm(const SpectralKSpace& spectral_kspace,
                                 const DistributionMapping& dm,
                                 const int norder_x, const int norder_y,
                                 const int norder_z, const Real dt);
        void pushSpectralFields( SpectralFieldData& f ) const;

    private:
        // Modified finite-order vectors
        SpectralKVector modified_kx_vec, modified_ky_vec, modified_kz_vec;
        SpectralCoefficients C_coef, S_ck_coef, X1_coef, X2_coef, X3_coef;
};

#endif // WARPX_PSATD_ALGORITHM_H_
