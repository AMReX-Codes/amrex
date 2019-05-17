
#include "MyTest.H"
#include <AMReX_PlotFileUtil.H>

using namespace amrex;

void
MyTest::writePlotfile () const
{
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
