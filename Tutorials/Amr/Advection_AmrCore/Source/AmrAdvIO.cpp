
#include <AMReX_PlotFileUtil.H>

#include <AmrAdv.H>

using namespace amrex;

std::string
AmrAdv::PlotFileName (int lev) const
{
    return amrex::Concatenate(plot_file, lev, 5);
}

Array<const MultiFab*>
AmrAdv::PlotFileMF () const
{
    Array<const MultiFab*> r;
    for (int i = 0; i <= finest_level; ++i) {
	r.push_back(phi_new[i].get());
    }
    return r;
}

Array<std::string>
AmrAdv::PlotFileVarNames () const
{
    return {"phi"};
}

void
AmrAdv::WritePlotFile () const
{
    const std::string& plotfilename = PlotFileName(istep[0]);
    const auto& mf = PlotFileMF();
    const auto& varnames = PlotFileVarNames();
    
    amrex::WriteMultiLevelPlotfile(plotfilename, finest_level+1, mf, varnames,
				    Geom(), t_new[0], istep, refRatio());
}
