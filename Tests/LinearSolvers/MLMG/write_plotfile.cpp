
#include <AMReX_PlotFileUtil.H>

using namespace amrex;

void write_plotfile (const Vector<Geometry>& geom, int rr,
                     const Vector<MultiFab>& soln, const Vector<MultiFab>& exact,
                     const Vector<MultiFab>& alpha, const Vector<MultiFab>& beta,
                     const Vector<MultiFab>& rhs)
{
    const int nlevels = geom.size();

    Vector<MultiFab> plotmf(nlevels);
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

        amrex::Print() << "Error on level " << ilev << ": " << plotmf[ilev].min(2)
                       << ", " << plotmf[ilev].max(2) << "\n";
    }

    amrex::WriteMultiLevelPlotfile ("plot",
                                    nlevels,
                                    amrex::GetVecOfConstPtrs(plotmf),
                                    {"solution", "exact", "error", "alpha", "beta", "rhs"},
                                    geom,
                                    0.0,
                                    Vector<int>(nlevels, 0),
                                    Vector<IntVect>(nlevels, IntVect{AMREX_D_DECL(rr,rr,rr)}));
}

