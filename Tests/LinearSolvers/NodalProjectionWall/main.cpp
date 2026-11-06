// Tests MLNodeLaplacian::setWallDivuCorrection, the second-order nodal divergence at walls
// and inflow faces, and shows the error it removes: after a nodal projection, a viscous
// operator applied to the velocity has an O(1) error in the wall-normal component in the cells
// next to a no-slip wall.  See README.md.
//
// checks = 1 (checks.cpp) checks compRHS without a solve, and stops with amrex::Abort if a
// check with a known answer fails.  study = 1 (study.cpp) projects test flows and prints
// the errors of the projected velocity and of the viscous operators applied to it.

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
        Real rtol = 1.e-12;
        pp.query("rtol", rtol);
        int verbose = 0;
        pp.query("verbose", verbose);
        int corr = 0;
        pp.query("wall_corr", corr);
        wall_corr = (corr != 0);
        amrex::Print() << "\nNodal projection wall correction (wall_corr): " << corr << "\n";

        int checks = 1;
        pp.query("checks", checks);
        int const nfail = (checks != 0) ? run_checks(n_cell, max_grid_size) : 0;

        int study = 0;
        pp.query("study", study);
        if (study != 0) {
            run_study(n_cell, max_grid_size, rtol, verbose);
        }

        if (nfail > 0) {
            amrex::Abort("NodalProjectionWall: " + std::to_string(nfail) + " checks failed");
        }
        amrex::Print() << "NodalProjectionWall: all checks passed\n";
    }
    amrex::Finalize();
}
