
#include <AMReX_PlotFileUtil.H>

using namespace amrex;

void write_plotfile (const Array<Geometry>& geom, int rr,
                     const Array<MultiFab>& soln, const Array<MultiFab>& exact,
                     const Array<MultiFab>& alpha, const Array<MultiFab>& beta,
                     const Array<MultiFab>& rhs)
{
    const int nlevels = geom.size();

    Array<MultiFab> plotmf(nlevels);
    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        plotmf[ilev].define(soln[ilev].boxArray(), soln[ilev].DistributionMap(), 6, 0);
        MultiFab::Copy(plotmf[ilev],  soln[ilev], 0, 0, 1, 0);
        MultiFab::Copy(plotmf[ilev], exact[ilev], 0, 1, 1, 0);
        MultiFab::Copy(plotmf[ilev],  soln[ilev], 0, 2, 1, 0);
        MultiFab::Copy(plotmf[ilev], alpha[ilev], 0, 3, 1, 0);
        MultiFab::Copy(plotmf[ilev],  beta[ilev], 0, 4, 1, 0);
        MultiFab::Copy(plotmf[ilev],   rhs[ilev], 0, 5, 1, 0);
        MultiFab::Subtract(plotmf[ilev], plotmf[ilev], 1, 2, 1, 0); // error = soln - exact
    }

    amrex::WriteMultiLevelPlotfile ("plot",
                                    nlevels,
                                    amrex::GetArrOfConstPtrs(plotmf),
                                    {"solution", "exact", "error", "alpha", "beta", "rhs"},
                                    geom,
                                    0.0,
                                    Array<int>(nlevels, 0),
                                    Array<IntVect>(nlevels, IntVect{AMREX_D_DECL(rr,rr,rr)}));
}

