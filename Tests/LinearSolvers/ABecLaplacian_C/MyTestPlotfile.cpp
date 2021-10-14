
#include "MyTest.H"
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

void
MyTest::writePlotfile () const
{
    ParmParse pp;
    bool gpu_regtest = false;
#ifdef AMREX_USE_GPU
    pp.query("gpu_regtest", gpu_regtest);
#endif

    const int nlevels = max_level+1;
    Vector<MultiFab> plotmf(nlevels);

    if (gpu_regtest) {
        const int ncomp = (acoef.empty()) ? 3 : 5;
        Vector<std::string> varname = {"solution", "rhs", "exact_solution"};
        if (!acoef.empty()) {
            varname.emplace_back("acoef");
            varname.emplace_back("bcoef");
        }

        for (int ilev = 0; ilev < nlevels; ++ilev)
        {
            plotmf[ilev].define(grids[ilev], dmap[ilev], ncomp, 0);
            MultiFab::Copy(plotmf[ilev], solution      [ilev], 0, 0, 1, 0);
            MultiFab::Copy(plotmf[ilev], rhs           [ilev], 0, 1, 1, 0);
            MultiFab::Copy(plotmf[ilev], exact_solution[ilev], 0, 2, 1, 0);
            if (!acoef.empty()) {
                MultiFab::Copy(plotmf[ilev], acoef[ilev], 0, 3, 1, 0);
                MultiFab::Copy(plotmf[ilev], bcoef[ilev], 0, 4, 1, 0);
            }
        }

        WriteMultiLevelPlotfile("plot", nlevels, amrex::GetVecOfConstPtrs(plotmf),
                                varname, geom, 0.0, Vector<int>(nlevels, 0),
                                Vector<IntVect>(nlevels, IntVect{ref_ratio}));
    } else {
        const int ncomp = (acoef.empty()) ? 4 : 6;
        Vector<std::string> varname = {"solution", "rhs", "exact_solution", "error"};
        if (!acoef.empty()) {
            varname.emplace_back("acoef");
            varname.emplace_back("bcoef");
        }

        for (int ilev = 0; ilev < nlevels; ++ilev)
        {
            plotmf[ilev].define(grids[ilev], dmap[ilev], ncomp, 0);
            MultiFab::Copy(plotmf[ilev], solution      [ilev], 0, 0, 1, 0);
            MultiFab::Copy(plotmf[ilev], rhs           [ilev], 0, 1, 1, 0);
            MultiFab::Copy(plotmf[ilev], exact_solution[ilev], 0, 2, 1, 0);
            MultiFab::Copy(plotmf[ilev], solution      [ilev], 0, 3, 1, 0);
            MultiFab::Subtract(plotmf[ilev], plotmf[ilev], 2, 3, 1, 0); // error = soln - exact
            if (!acoef.empty()) {
                MultiFab::Copy(plotmf[ilev], acoef[ilev], 0, 4, 1, 0);
                MultiFab::Copy(plotmf[ilev], bcoef[ilev], 0, 5, 1, 0);
            }
            auto dx = geom[ilev].CellSize();
            Real dvol = AMREX_D_TERM(dx[0],*dx[1],*dx[2]);
            amrex::Print() << "Level " << ilev
                           << " max-norm error: " << plotmf[ilev].norminf(3)
                           << " 1-norm error: " << plotmf[ilev].norm1(3)*dvol << std::endl;
        }

        WriteMultiLevelPlotfile("plot", nlevels, amrex::GetVecOfConstPtrs(plotmf),
                                varname, geom, 0.0, Vector<int>(nlevels, 0),
                                Vector<IntVect>(nlevels, IntVect{ref_ratio}));
    }
}

