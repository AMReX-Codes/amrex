
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
    if (gpu_regtest) {
        const int ncomp = AMREX_SPACEDIM*3 + 2;
        Vector<std::string> varname =
            {"u", "v", "w", "uexact", "vexact", "wexact",
             "xrhs", "yrhs", "zrhs", "eta", "vfrc"};

        const MultiFab& vfrc = factory->getVolFrac();

        MultiFab plotmf(grids, dmap, ncomp,  0);
        MultiFab::Copy(plotmf, solution, 0,  0, 3, 0);
        MultiFab::Copy(plotmf, exact   , 0,  3, 3, 0);
        MultiFab::Copy(plotmf, rhs     , 0,  6, 3, 0);
        MultiFab::Copy(plotmf, eta     , 0,  9, 1, 0);
        MultiFab::Copy(plotmf, vfrc    , 0, 10, 1, 0);

        WriteMultiLevelPlotfile("plot", 1, {&plotmf},
                                varname, {geom}, 0.0, {0}, {IntVect(2)});
    } else {
        const int ncomp = AMREX_SPACEDIM*4 + 2;
        Vector<std::string> varname =
            {"u", "v", "w", "uexact", "vexact", "wexact", "xerror", "yerror", "zerror",
             "xrhs", "yrhs", "zrhs", "eta", "vfrc"};

        const MultiFab& vfrc = factory->getVolFrac();

        MultiFab plotmf(grids, dmap, ncomp,  0);
        MultiFab::Copy(plotmf, solution, 0,  0, 3, 0);
        MultiFab::Copy(plotmf, exact   , 0,  3, 3, 0);
        MultiFab::Copy(plotmf, solution, 0,  6, 3, 0);
        MultiFab::Copy(plotmf, rhs     , 0,  9, 3, 0);
        MultiFab::Copy(plotmf, eta     , 0, 12, 1, 0);
        MultiFab::Copy(plotmf, vfrc    , 0, 13, 1, 0);

        MultiFab::Subtract(plotmf, exact, 0, 6, 3, 0);

        WriteMultiLevelPlotfile("plot", 1, {&plotmf},
                                varname, {geom}, 0.0, {0}, {IntVect(2)});
    }
}
