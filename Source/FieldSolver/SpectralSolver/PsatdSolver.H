#ifndef WARPX_PSATD_SOLVER_H_
#define WARPX_PSATD_SOLVER_H_

#include <AMREX_Real.H>
#include <WarpX_Complex.H>
#include <SpectralData.H>

using namespace amrex;
using namespace Gpu;

using SpectralVector = LayoutData<ManagedVector<Real>>;

/* TODO: Write documentation
*/
class PsatdSolver
{
    using SpectralCoefficients = FabArray<BaseFab<Real>>;

    public:
        PsatdSolver( const BoxArray& ba, const DistributionMapping& dm,
                     const Real* dx, const Real dt );
        void pushSpectralFields( SpectralData& f ) const;

    private:
        SpectralVector kx_vec, ky_vec, kz_vec;
        SpectralCoefficients C_coef, S_ck_coef, X1_coef, X2_coef, X3_coef;
};

#endif // WARPX_PSATD_SOLVER_H_
