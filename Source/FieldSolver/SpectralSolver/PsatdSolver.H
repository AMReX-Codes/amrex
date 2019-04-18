#ifndef WARPX_PSATD_SOLVER_H_
#define WARPX_PSATD_SOLVER_H_

#include <AMREX_Real.H>
#include <WarpXComplex.H>
#include <SpectralData.H>

using namespace amrex;
using namespace Gpu;

using SpectralVector = LayoutData<ManagedVector<Real>>

class PsatdSolver
{
    using SpectralCoefficients = FabArray<BaseFab<Real>>

    public:
        PsatdSolver( const BoxArray& ba, const DistributionMapping& dm, const Real* dx );
        void pushSpectralFields( SpectralData& f ) const;

    private:
        SpectralVector kx, ky, kz;
        SpectralCoefficients C_coef, S_ck_coef, X1_coef, X2_coef, X3_coef;
};

#endif // WARPX_PSATD_SOLVER_H_
