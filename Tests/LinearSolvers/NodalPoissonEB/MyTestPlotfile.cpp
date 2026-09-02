
#include "MyTest.H"
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>

using namespace amrex;

void
MyTest::writePlotfile () const
{
    Vector<std::string> varname;
    if (gpu_regtest) {
        varname = Vector<std::string>{"solution", "rhs", "exact_solution"};
    } else {
        varname = Vector<std::string>{"solution", "rhs", "exact_solution", "error"};
    }
    const auto ncomp = int(varname.size());

    MultiFab plotmf(grids[0], dmap[0], ncomp, 0, MFInfo(), *factory);
    amrex::average_node_to_cellcenter(plotmf, 0, solution[0], 0, 1);
    amrex::average_node_to_cellcenter(plotmf, 1, rhs[0], 0, 1);
    amrex::average_node_to_cellcenter(plotmf, 2, exact_solution[0], 0, 1);

    if (!gpu_regtest) {
        MultiFab error(rhs[0].boxArray(), rhs[0].DistributionMap(), 1, 0);
        MultiFab::Copy(error, solution[0], 0, 0, 1, 0);
        MultiFab::Subtract(error, exact_solution[0], 0, 0, 1, 0);
        amrex::average_node_to_cellcenter(plotmf, 3, error, 0, 1);
    }

    EB_WriteSingleLevelPlotfile("plot", plotmf, varname, geom[0], 0.0, 0);
}
