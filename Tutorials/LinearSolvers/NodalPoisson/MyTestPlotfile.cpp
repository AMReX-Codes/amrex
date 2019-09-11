
#include "MyTest.H"
#include <AMReX_PlotFileUtil.H>

using namespace amrex;

void
MyTest::writePlotfile () const
{
    const int ncomp = 4;
    Vector<std::string> varname = {"solution", "rhs", "exact_solution", "error"};

    const int nlevels = max_level+1;

    Vector<MultiFab> plotmf(nlevels);
    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        plotmf[ilev].define(grids[ilev], dmap[ilev], ncomp, 0);
        MultiFab::Copy(plotmf[ilev], solution      [ilev], 0, 0, 1, 0);
        MultiFab::Copy(plotmf[ilev], rhs           [ilev], 0, 1, 1, 0);
        MultiFab::Copy(plotmf[ilev], exact_solution[ilev], 0, 2, 1, 0);
        MultiFab::Copy(plotmf[ilev], solution      [ilev], 0, 3, 1, 0);
        MultiFab::Subtract(plotmf[ilev], plotmf[ilev], 2, 3, 1, 0); // error = soln - exact
    }

    WriteMultiLevelPlotfile("plot", nlevels, amrex::GetVecOfConstPtrs(plotmf),
                            varname, geom, 0.0, Vector<int>(nlevels, 0),
                            Vector<IntVect>(nlevels, IntVect{ref_ratio}));
}

