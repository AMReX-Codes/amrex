
#include <AMReX_MultiFab.H>
#include <AMReX_MLMG.H>

using namespace amrex;

void solve_with_mlmg (const Array<Geometry>& geom, Array<MultiFab>& soln,
                      const Array<MultiFab>& alpha, const Array<MultiFab>& beta,
                      const Array<MultiFab>& rhs)
{
    const Real tol_rel = 1.e-10;
    const Real tol_abs = 0.0;

}

