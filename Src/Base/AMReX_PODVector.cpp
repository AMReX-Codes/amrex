#include <AMReX_PODVector.H>
#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>

namespace amrex::VectorGrowthStrategy
{
    Real growth_factor = 1.5._rt;

    void Initialize () {
        ParmParse pp("amrex");
        pp.queryAdd("vector_growth_factor", growth_factor);
    }
}
