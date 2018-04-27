#include <AMReX.H>
#include <AMReX_PlotFileUtil.H>

using namespace amrex;

extern "C"
{
    void amrex_fi_write_plotfile (const char* name, int nlevs, const MultiFab* mf[],
				  const char* varname[], const Geometry* geom[],
				  Real time, const int level_steps[], const int ref_ratio[])
    {
	Vector<const MultiFab*> mfarr {mf, mf+nlevs};
	Vector<std::string> varnamearr {varname, varname+mf[0]->nComp()};
	Vector<Geometry> geomarr;
	for (int lev = 0; lev < nlevs; ++lev) {
	    geomarr.emplace_back(*geom[lev]);
	}
	Vector<int> lsarr {level_steps, level_steps+nlevs};
	Vector<IntVect> rrarr;
	for (int lev = 0; lev < nlevs-1; ++lev) {
	    rrarr.emplace_back(AMREX_D_DECL(ref_ratio[lev],ref_ratio[lev],ref_ratio[lev]));
	}
	amrex::WriteMultiLevelPlotfile(name, nlevs, mfarr, varnamearr, geomarr,
				       time, lsarr, rrarr);
    }
}
