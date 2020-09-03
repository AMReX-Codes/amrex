
#include "MyTest.H"
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

void
MyTest::writePlotfile () const
{
    ParmParse pp;
    bool plot_error = true;

    pp.query("plot_error", plot_error);

    if (plot_error) {
        const int ncomp = AMREX_SPACEDIM*4 + 2;
        Vector<std::string> varname =
#if (AMREX_SPACEDIM == 2)
            {"u", "v", "uexact", "vexact", "xerror", "yerror", 
             "xrhs", "yrhs", "eta", "vfrc"};
#else
            {"u", "v", "w", "uexact", "vexact", "wexact", "xerror", "yerror", "zerror",
             "xrhs", "yrhs", "zrhs", "eta", "vfrc"};
#endif

        const MultiFab& vfrc = factory->getVolFrac();

        MultiFab plotmf(grids, dmap, ncomp,  0);
        MultiFab::Copy(plotmf, solution, 0,                 0  , AMREX_SPACEDIM, 0);
        MultiFab::Copy(plotmf, exact   , 0,    AMREX_SPACEDIM  , AMREX_SPACEDIM, 0);
        MultiFab::Copy(plotmf, solution, 0,  2*AMREX_SPACEDIM  , AMREX_SPACEDIM, 0);
        MultiFab::Copy(plotmf, rhs     , 0,  3*AMREX_SPACEDIM  , AMREX_SPACEDIM, 0);
        MultiFab::Copy(plotmf, eta     , 0,  4*AMREX_SPACEDIM  , 1, 0);
        MultiFab::Copy(plotmf, vfrc    , 0,  4*AMREX_SPACEDIM+1, 1, 0);

        MultiFab::Subtract(plotmf, exact, 0, 2*AMREX_SPACEDIM, AMREX_SPACEDIM, 0);

        WriteMultiLevelPlotfile("plot", 1, {&plotmf},
                                varname, {geom}, 0.0, {0}, {IntVect(2)});

    } else {

        const int ncomp = AMREX_SPACEDIM*3 + 2;
        Vector<std::string> varname =
#if (AMREX_SPACEDIM == 2)
            {"u", "v", "uexact", "vexact",
             "xrhs", "yrhs", "eta", "vfrc"};
#else
            {"u", "v", "w", "uexact", "vexact", "wexact",
             "xrhs", "yrhs", "zrhs", "eta", "vfrc"};
#endif

        const MultiFab& vfrc = factory->getVolFrac();

        MultiFab plotmf(grids, dmap, ncomp,  0);
        MultiFab::Copy(plotmf, solution, 0,               0  , AMREX_SPACEDIM, 0);
        MultiFab::Copy(plotmf, exact   , 0,  AMREX_SPACEDIM  , AMREX_SPACEDIM, 0);
        MultiFab::Copy(plotmf, rhs     , 0,2*AMREX_SPACEDIM  , AMREX_SPACEDIM, 0);
        MultiFab::Copy(plotmf, eta     , 0,3*AMREX_SPACEDIM  , 1, 0);
        MultiFab::Copy(plotmf, vfrc    , 0,3*AMREX_SPACEDIM+1, 1, 0);

        WriteMultiLevelPlotfile("plot", 1, {&plotmf},
                                varname, {geom}, 0.0, {0}, {IntVect(2)});
    }
}
