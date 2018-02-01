#include "MyTest.H"

#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_EBIndexSpace.H>
#include <AMReX_EBTower.H>

using namespace amrex;

MyTest::MyTest ()
{
    readParameters();

    initGrids();

    initializeEBIS();

    EBTower::Build();
    AMReX_EBIS::reset();

    initData();
}

//
// Given vel, rhs & sig, this solves Div (sig * Grad phi) = Div vel + rhs.
// On return, vel becomes vel  - sig * Grad phi.
//
void
MyTest::solve ()
{
    std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_lobc;
    std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_hibc;
    mlmg_lobc[0] = LinOpBCType::Neumann;
    mlmg_hibc[0] = LinOpBCType::Neumann;
    mlmg_lobc[1] = LinOpBCType::Neumann;
    mlmg_hibc[1] = LinOpBCType::Dirichlet;
    static_assert(AMREX_SPACEDIM==2, "2d only");

    MLNodeLaplacian mlndlap(geom, grids, dmap);

    mlndlap.setDomainBC(mlmg_lobc, mlmg_hibc);

    for (int ilev = 0; ilev <= max_level; ++ilev) {
        mlndlap.setSigma(ilev, sig[ilev]);
    }

    mlndlap.compRHS(amrex::GetVecOfPtrs(rhs), amrex::GetVecOfPtrs(vel), {}, {});

    MLMG mlmg(mlndlap);
    mlmg.setVerbose(verbose);
    mlmg.setBottomVerbose(bottom_verbose);
    mlmg.setMaxIter(max_iter);
    mlmg.setMaxFmgIter(max_fmg_iter);

    Real mlmg_err = mlmg.solve(amrex::GetVecOfPtrs(phi), amrex::GetVecOfConstPtrs(rhs),
                               1.e-11, 0.0);

    mlndlap.updateVelocity(amrex::GetVecOfPtrs(vel), amrex::GetVecOfConstPtrs(phi));

    for (int ilev = 0; ilev <= max_level; ++ilev) {
        amrex::VisMF::Write(phi[ilev], "phi-"+std::to_string(ilev));
        amrex::VisMF::Write(rhs[ilev], "rhs-"+std::to_string(ilev));
        amrex::VisMF::Write(vel[ilev], "vel-"+std::to_string(ilev));
        amrex::VisMF::Write(factory[ilev]->getVolFrac(), "vfrc-"+std::to_string(ilev));
    }
}

void
MyTest::readParameters ()
{
    ParmParse pp;
    pp.query("max_level", max_level);
    pp.query("n_cell", n_cell);
    pp.query("max_grid_size", max_grid_size);

    pp.query("verbose", verbose);
    pp.query("bottom_verbose", bottom_verbose);
    pp.query("max_iter", max_iter);
    pp.query("max_fmg_iter", max_fmg_iter);
#ifdef AMREX_USE_HYPRE
    pp.query("use_hypre", use_hypre);
#endif

    pp.query("geom_type", geom_type);
}

void
MyTest::initGrids ()
{
    int nlevels = max_level + 1;
    geom.resize(nlevels);
    grids.resize(nlevels);

    RealBox rb({AMREX_D_DECL(0.,0.,0.)}, {AMREX_D_DECL(1.,1.,1.)});
    std::array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0,0,0)};
    Geometry::Setup(&rb, 0, is_periodic.data());
    Box domain0(IntVect{AMREX_D_DECL(0,0,0)}, IntVect{AMREX_D_DECL(n_cell-1,n_cell-1,n_cell-1)});
    Box domain = domain0;
    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        geom[ilev].define(domain);
        domain.refine(ref_ratio);
    }

    domain = domain0;
    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        grids[ilev].define(domain);
        grids[ilev].maxSize(max_grid_size);
        domain.grow(-n_cell/4);   // fine level cover the middle of the coarse domain
        domain.refine(ref_ratio); 
    }
}

void
MyTest::initData ()
{
    int nlevels = max_level + 1;
    dmap.resize(nlevels);
    factory.resize(nlevels);
    phi.resize(nlevels);
    rhs.resize(nlevels);
    vel.resize(nlevels);
    sig.resize(nlevels);

    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        dmap[ilev].define(grids[ilev]);
        factory[ilev].reset(new EBFArrayBoxFactory(geom[ilev], grids[ilev], dmap[ilev],
                                                   {2,2,2}, EBSupport::volume));

        phi[ilev].define(amrex::convert(grids[ilev],IntVect::TheNodeVector()),
                         dmap[ilev], 1, 1, MFInfo(), *factory[ilev]);
        rhs[ilev].define(amrex::convert(grids[ilev],IntVect::TheNodeVector()),
                         dmap[ilev], 1, 0, MFInfo(), *factory[ilev]);
        sig[ilev].define(grids[ilev], dmap[ilev], 1, 1, MFInfo(), *factory[ilev]);
        vel[ilev].define(grids[ilev], dmap[ilev], AMREX_SPACEDIM, 1, MFInfo(), *factory[ilev]);

        phi[ilev].setVal(0.0);
        sig[ilev].setVal(1.0);
        vel[ilev].setVal(0.0);
        const int icomp = 1; // vy
        const int ncomp = 1;
        for (MFIter mfi(vel[ilev]); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.validbox();
            FArrayBox& fab = vel[ilev][mfi];
            fab.ForEachIV(bx, icomp, ncomp, [=] (Real& x, const IntVect& iv)
                          { if (iv[1] < n_cell/2) x = 1.0; });
        }
    }

}

