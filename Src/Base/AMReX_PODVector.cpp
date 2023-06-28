#include <AMReX_PODVector.H>
#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>

namespace amrex::VectorGrowthStrategy
{
    Real growth_factor = 1.5_rt;

    void Initialize () {
        ParmParse pp("amrex");
        pp.queryAdd("vector_growth_factor", growth_factor);

        // sanity checks
        auto eps = std::numeric_limits<Real>::epsilon();
        auto huge = 1000._rt;  // huge enough for our purposes...
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(growth_factor - 1.0_rt >= eps,
            "User-specified vector growth factor is too small.");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(growth_factor < huge,
            "User-specified vector growth factor is too large.");
    }
}
