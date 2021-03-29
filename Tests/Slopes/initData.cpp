#include "MyTest.H"
#include <AMReX_Config.H>
#include <AMReX_EB2.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_MultiFabUtil.H>

using namespace amrex;

void
MyTest::initData ()
{
    int nlevels = max_level + 1;
    dmap.resize(nlevels);
    factory.resize(nlevels);
    phi.resize(nlevels);
    grad_x.resize(nlevels);
    grad_x_analytic.resize(nlevels);
    grad_y.resize(nlevels);
    grad_y_analytic.resize(nlevels);
    grad_z.resize(nlevels);
    grad_z_analytic.resize(nlevels);
    ccentr.resize(nlevels);

    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        dmap[ilev].define(grids[ilev]);
        const EB2::IndexSpace& eb_is = EB2::IndexSpace::top();
        const EB2::Level& eb_level = eb_is.getLevel(geom[ilev]);
        factory[ilev].reset(new EBFArrayBoxFactory(eb_level, geom[ilev], grids[ilev], dmap[ilev],
                                                   {2,2,2}, EBSupport::full));

        phi[ilev].define(             grids[ilev], dmap[ilev], AMREX_SPACEDIM, 1, MFInfo(), *factory[ilev]);
        grad_x[ilev].define(          grids[ilev], dmap[ilev], AMREX_SPACEDIM, 1, MFInfo(), *factory[ilev]);
        grad_x_analytic[ilev].define( grids[ilev], dmap[ilev], AMREX_SPACEDIM, 1, MFInfo(), *factory[ilev]);
        grad_y[ilev].define(          grids[ilev], dmap[ilev], AMREX_SPACEDIM, 1, MFInfo(), *factory[ilev]);
        grad_y_analytic[ilev].define( grids[ilev], dmap[ilev], AMREX_SPACEDIM, 1, MFInfo(), *factory[ilev]);
        grad_z[ilev].define(          grids[ilev], dmap[ilev], AMREX_SPACEDIM, 1, MFInfo(), *factory[ilev]);
        grad_z_analytic[ilev].define( grids[ilev], dmap[ilev], AMREX_SPACEDIM, 1, MFInfo(), *factory[ilev]);
        ccentr[ilev].define(          grids[ilev], dmap[ilev], AMREX_SPACEDIM, 0, MFInfo(), *factory[ilev]);

        phi[ilev].setVal(0.0);
        grad_x[ilev].setVal(1e40);
        grad_x_analytic[ilev].setVal(1e40);
        grad_y[ilev].setVal(1e40);
        grad_y_analytic[ilev].setVal(1e40);
        grad_z[ilev].setVal(1e40);
        grad_z_analytic[ilev].setVal(1e40);
        ccentr[ilev].setVal(0.0);

        if(use_linear_1d)
        {
            initializeLinearData(ilev);

        } else {
            for (MFIter mfi(phi[ilev]); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.fabbox();
                Array4<Real> const& fab = phi[ilev].array(mfi);
                const auto dx = geom[ilev].CellSizeArray();

                amrex::ParallelFor(bx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    Real rx = (i+0.5)*dx[0];
                    Real ry = (j+0.5)*dx[1];
                    fab(i,j,k) = rx*(1.-rx)*ry*(1.-ry);
                });
            }
        }
    }
}
