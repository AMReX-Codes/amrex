#include "MyTest.H"
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
    phi_ghost_resolved.resize(nlevels);
    phieb.resize(nlevels);
    grad_x.resize(nlevels);
    grad_x_analytic.resize(nlevels);
    grad_y.resize(nlevels);
    grad_y_analytic.resize(nlevels);
    grad_z.resize(nlevels);
    grad_z_analytic.resize(nlevels);
    grad_eb.resize(nlevels);
    grad_eb_analytic.resize(nlevels);
    lap_analytic.resize(nlevels);
    ccentr.resize(nlevels);
    rhs.resize(nlevels);
    acoef.resize(nlevels);
    bcoef.resize(nlevels);
    bcoef_eb.resize(nlevels);


    m_bcrec.resize(1);
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      if (geom[0].isPeriodic(dir)) {
        m_bcrec[0].setLo(dir, BCType::int_dir);
        m_bcrec[0].setHi(dir, BCType::int_dir);
      } else {
        m_bcrec[0].setLo(dir, BCType::ext_dir);
        m_bcrec[0].setHi(dir, BCType::ext_dir);
      }
    }


    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        dmap[ilev].define(grids[ilev]);
        const EB2::IndexSpace& eb_is = EB2::IndexSpace::top();
        const EB2::Level& eb_level = eb_is.getLevel(geom[ilev]);
        factory[ilev] = std::make_unique<EBFArrayBoxFactory>
            (eb_level, geom[ilev], grids[ilev], dmap[ilev], Vector<int>{2,2,2}, EBSupport::full);

        phi[ilev].define(grids[ilev], dmap[ilev], AMREX_SPACEDIM, 1, MFInfo(), *factory[ilev]);
        phi_ghost_resolved[ilev].define(grids[ilev], dmap[ilev], AMREX_SPACEDIM, 1, MFInfo(), *factory[ilev]);
        phieb[ilev].define(grids[ilev], dmap[ilev], AMREX_SPACEDIM, 1, MFInfo(), *factory[ilev]);
        grad_x[ilev].define(grids[ilev], dmap[ilev], AMREX_SPACEDIM, 1, MFInfo(), *factory[ilev]);
        grad_x_analytic[ilev].define(grids[ilev], dmap[ilev], AMREX_SPACEDIM, 1, MFInfo(), *factory[ilev]);
        grad_y[ilev].define(grids[ilev], dmap[ilev], AMREX_SPACEDIM, 1, MFInfo(), *factory[ilev]);
        grad_y_analytic[ilev].define(grids[ilev], dmap[ilev], AMREX_SPACEDIM, 1, MFInfo(), *factory[ilev]);
        grad_z[ilev].define(grids[ilev], dmap[ilev], AMREX_SPACEDIM, 1, MFInfo(), *factory[ilev]);
        grad_z_analytic[ilev].define(grids[ilev], dmap[ilev], AMREX_SPACEDIM, 1, MFInfo(), *factory[ilev]);
        grad_eb[ilev].define(grids[ilev], dmap[ilev], AMREX_SPACEDIM, 1, MFInfo(), *factory[ilev]);
        grad_eb_analytic[ilev].define(grids[ilev], dmap[ilev], AMREX_SPACEDIM, 1, MFInfo(), *factory[ilev]);
        lap_analytic[ilev].define(grids[ilev], dmap[ilev], AMREX_SPACEDIM, 1, MFInfo(), *factory[ilev]);
        ccentr[ilev].define(grids[ilev], dmap[ilev], AMREX_SPACEDIM, 1, MFInfo(), *factory[ilev]);
        rhs[ilev].define(grids[ilev], dmap[ilev], AMREX_SPACEDIM, 0, MFInfo(), *factory[ilev]);
        acoef[ilev].define(grids[ilev], dmap[ilev], AMREX_SPACEDIM, 0, MFInfo(), *factory[ilev]);
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            bcoef[ilev][idim].define(amrex::convert(grids[ilev],IntVect::TheDimensionVector(idim)),
                                     dmap[ilev], AMREX_SPACEDIM, 0, MFInfo(), *factory[ilev]);
        }
        if (eb_is_dirichlet || eb_is_homog_dirichlet) {
            bcoef_eb[ilev].define(grids[ilev], dmap[ilev], AMREX_SPACEDIM, 0, MFInfo(), *factory[ilev]);
            bcoef_eb[ilev].setVal(1.0);
        }

        phi[ilev].setVal(0.0);
        phi_ghost_resolved[ilev].setVal(0.0);
        phieb[ilev].setVal(0.0);
        grad_x[ilev].setVal(1e40);
        grad_x_analytic[ilev].setVal(1e40);
        grad_y[ilev].setVal(1e40);
        grad_y_analytic[ilev].setVal(1e40);
        grad_z[ilev].setVal(1e40);
        grad_z_analytic[ilev].setVal(1e40);
        grad_eb[ilev].setVal(1e40);
        grad_eb_analytic[ilev].setVal(1e40);
        lap_analytic[ilev].setVal(1e40);
        ccentr[ilev].setVal(0.0);
        rhs[ilev].setVal(0.0);
        acoef[ilev].setVal(0.0);
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            bcoef[ilev][idim].setVal(1.0);
        }

        if(use_poiseuille) {
           initializePoiseuilleData(ilev);
        }
        else if(use_triangle_wave) {
           initializeTriangleWaveData(ilev);
        }
        else {
            // Test a custom polynomial function
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
                    // fab(i,j,k) = std::sin(rx*2.*pi + 43.5)*std::sin(ry*2.*pi + 89.);
                    fab(i,j,k) = rx*(1.-rx)*ry*(1.-ry);
                });
            }
        }
    }
}
