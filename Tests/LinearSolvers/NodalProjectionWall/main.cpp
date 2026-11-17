// Tests MLNodeLaplacian::setBoundaryDivuCorrection, the second-order nodal divergence at walls
// and inflow faces: checks of compRHS without a solve (checks.cpp), which stop the test with
// amrex::Abort if a check with a known answer fails.  See README.md.

#include "NodalProjectionWall.H"

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <string>

using namespace amrex;

bool wall_corr = false;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        ParmParse pp;
        // queryarr does not shrink the vector, so the defaults are set only if n_cell is absent
        Vector<int> n_cell;
        if (!pp.queryarr("n_cell", n_cell)) {
            n_cell = {32, 64, 128, 256};
        }
        int max_grid_size = 32;
        pp.query("max_grid_size", max_grid_size);
        int corr = 0;
        pp.query("wall_corr", corr);
        wall_corr = (corr != 0);
        amrex::Print() << "\nNodal projection wall correction (wall_corr): " << corr << "\n";

        int const nfail = run_checks(n_cell, max_grid_size);

        if (nfail > 0) {
            amrex::Abort("NodalProjectionWall: " + std::to_string(nfail) + " checks failed");
        }
        amrex::Print() << "NodalProjectionWall: all checks passed\n";
    }
    amrex::Finalize();
}
