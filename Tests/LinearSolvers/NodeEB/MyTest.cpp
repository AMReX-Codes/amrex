#include "MyTest.H"

#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_EB2.H>

#if(AMREX_SPACEDIM == 3)
#include <AMReX_algoim_integrals.H>
#endif

using namespace amrex;

MyTest::MyTest ()
{
    readParameters();

    initGrids();

    initializeEB();

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
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (Geometry::isPeriodic(idim)) {
            mlmg_lobc[idim] = LinOpBCType::Periodic;
            mlmg_hibc[idim] = LinOpBCType::Periodic;
        } else {
            mlmg_lobc[idim] = LinOpBCType::Neumann;
            mlmg_hibc[idim] = LinOpBCType::Dirichlet;
        }
    }
            

    LPInfo info;
    info.setMaxCoarseningLevel(max_coarsening_level);

    MLNodeLaplacian mlndlap(geom, grids, dmap, info, amrex::GetVecOfConstPtrs(factory));

    if (sigma) {
        mlndlap.setCoarseningStrategy(MLNodeLaplacian::CoarseningStrategy::Sigma);
    }

    mlndlap.setDomainBC(mlmg_lobc, mlmg_hibc);

    for (int ilev = 0; ilev <= max_level; ++ilev) {
        mlndlap.setSigma(ilev, sig[ilev]);
    }

    mlndlap.compRHS(amrex::GetVecOfPtrs(rhs), amrex::GetVecOfPtrs(vel), {}, {});

#if 0
#if (AMREX_SPACEDIM == 2)
    for (int ilev = 0; ilev <= max_level; ++ilev) {
        amrex::VisMF::Write(rhs[ilev], "rhs2d");
    }
#else
    for (int ilev = 0; ilev <= max_level; ++ilev) {
        amrex::VisMF::Write(rhs[ilev], "rhs3d");
    }
#endif
#endif

    MLMG mlmg(mlndlap);
    mlmg.setVerbose(verbose);
    mlmg.setBottomVerbose(bottom_verbose);
    mlmg.setMaxIter(max_iter);
    mlmg.setMaxFmgIter(max_fmg_iter);

#ifdef AMREX_USE_HYPRE
    if (use_hypre) mlmg.setBottomSolver(BottomSolver::hypre);
#endif

    Real mlmg_err = mlmg.solve(amrex::GetVecOfPtrs(phi), amrex::GetVecOfConstPtrs(rhs),
                               1.e-11, 0.0);

    mlndlap.updateVelocity(amrex::GetVecOfPtrs(vel), amrex::GetVecOfConstPtrs(phi));

#if 0
#if (AMREX_SPACEDIM == 2)
    for (int ilev = 0; ilev <= max_level; ++ilev) {
        amrex::VisMF::Write(phi[ilev], "phi2d");
    }
#else
    for (int ilev = 0; ilev <= max_level; ++ilev) {
        amrex::VisMF::Write(phi[ilev], "phi3d");
    }
#endif
#endif

    mlndlap.compRHS(amrex::GetVecOfPtrs(rhs), amrex::GetVecOfPtrs(vel), {}, {});

#if 0
    for (int ilev = 0; ilev <= max_level; ++ilev) {
        amrex::VisMF::Write(rhs[ilev], "rhs"+std::to_string(ilev));
        amrex::Print() << "rhs.norm0() = " << rhs[ilev].norm0() << "\n";
        amrex::Print() << "rhs.norm1()/npoints = " << rhs[ilev].norm1() / grids[0].d_numPts() << "\n";
    }
#endif
}

void
MyTest::writePlotfile ()
{
    amrex::WriteSingleLevelPlotfile("plot", vel[0], {AMREX_D_DECL("xvel","yvel","zvel")}, geom[0], 0.0, 0);
}

void
MyTest::readParameters ()
{   
    ParmParse pp;
    pp.query("max_level", max_level);
    pp.query("n_cell", n_cell);
    pp.query("max_grid_size", max_grid_size);

    pp.query("plot_file", plot_file_name);

    pp.query("verbose", verbose);
    pp.query("bottom_verbose", bottom_verbose);
    pp.query("max_iter", max_iter);
    pp.query("max_fmg_iter", max_fmg_iter);
    pp.query("max_coarsening_level", max_coarsening_level);
#ifdef AMREX_USE_HYPRE
    pp.query("use_hypre", use_hypre);
#endif

    pp.query("sigma", sigma);

#if (AMREX_SPACEDIM == 3)
    ParmParse pp_eb("eb2");
    std::string geom_type;
    pp_eb.get("cylinder_direction", cylinder_direction);
#endif
}

void
MyTest::initGrids ()
{
    int nlevels = max_level + 1;
    geom.resize(nlevels);
    grids.resize(nlevels);

    RealBox rb({AMREX_D_DECL(0.,0.,0.)}, {AMREX_D_DECL(1.,1.,1.)});

    // Make the domain periodic at the ends of the cylinder
    if (cylinder_direction == 0)
    {
       std::array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1,0,0)};
       Geometry::Setup(&rb, 0, is_periodic.data());

    } else if (cylinder_direction == 1)
    {
       std::array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0,1,0)};
       Geometry::Setup(&rb, 0, is_periodic.data());

    } else if (cylinder_direction == 2)
    {
       std::array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0,0,1)};
       Geometry::Setup(&rb, 0, is_periodic.data());
    }

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
        const EB2::IndexSpace& eb_is = EB2::IndexSpace::top();
        const EB2::Level& eb_level = eb_is.getLevel(geom[ilev]);
        factory[ilev].reset(new EBFArrayBoxFactory(eb_level, geom[ilev], grids[ilev], dmap[ilev],
                                                   {2,2,2}, EBSupport::full));

        phi[ilev].define(amrex::convert(grids[ilev],IntVect::TheNodeVector()),
                         dmap[ilev], 1, 1, MFInfo(), *factory[ilev]);
        rhs[ilev].define(amrex::convert(grids[ilev],IntVect::TheNodeVector()),
                         dmap[ilev], 1, 0, MFInfo(), *factory[ilev]);
        sig[ilev].define(grids[ilev], dmap[ilev], 1, 1, MFInfo(), *factory[ilev]);
        vel[ilev].define(grids[ilev], dmap[ilev], AMREX_SPACEDIM, 1, MFInfo(), *factory[ilev]);

        phi[ilev].setVal(0.0);
        sig[ilev].setVal(1.0);
        vel[ilev].setVal(0.0);
        const Real* dx = geom[ilev].CellSize();
        const Real h = dx[0];
        for (MFIter mfi(vel[ilev]); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.validbox();
            FArrayBox& fab = vel[ilev][mfi];

            int ncomp = 1;

#if (AMREX_SPACEDIM > 2)
            if (cylinder_direction == 2)
#endif
            {
               int icomp = 0; // vx
               fab.ForEachIV(bx, icomp, ncomp, [=] (Real& v_x, const IntVect& iv) 
                   {
                       Real rx = (iv[0]+0.5)*h - 0.5;
                       Real ry = (iv[1]+0.5)*h - 0.5;
                       Real r = std::sqrt(rx*rx+ry*ry);
                       Real fac = std::exp(-(r*r/(0.16*0.16)));
                       v_x = v_x + 2.0*r*ry/r*fac;
                   });
               icomp = 1; // vy
               fab.ForEachIV(bx, icomp, ncomp, [=] (Real& v_y, const IntVect& iv) 
                   {
                       Real rx = (iv[0]+0.5)*h - 0.5;
                       Real ry = (iv[1]+0.5)*h - 0.5;
                       Real r = std::sqrt(rx*rx+ry*ry);
                       Real fac = std::exp(-(r*r/(0.16*0.16)));
                       v_y = v_y - 2.0*r*rx/r*fac;
                   });

            } 
#if (AMREX_SPACEDIM > 2)
            else if (cylinder_direction == 1) 
            {
               int icomp = 0; // vx
               fab.ForEachIV(bx, icomp, ncomp, [=] (Real& v_x, const IntVect& iv) 
                   {
                       Real rx = (iv[0]+0.5)*h - 0.5;
                       Real rz = (iv[2]+0.5)*h - 0.5;
                       Real r = std::sqrt(rx*rx+rz*rz);
                       Real fac = std::exp(-(r*r/(0.16*0.16)));
                       v_x = v_x - 2.0*r*rz/r*fac;
                   });
               icomp = 2; // vx
               fab.ForEachIV(bx, icomp, ncomp, [=] (Real& v_z, const IntVect& iv) 
                   {
                       Real rx = (iv[0]+0.5)*h - 0.5;
                       Real rz = (iv[2]+0.5)*h - 0.5;
                       Real r = std::sqrt(rx*rx+rz*rz);
                       Real fac = std::exp(-(r*r/(0.16*0.16)));
                       v_z = v_z + 2.0*r*rx/r*fac;
                   });

            } 
            else if (cylinder_direction == 0) 
            {
               int icomp = 1; // v:
               fab.ForEachIV(bx, icomp, ncomp, [=] (Real& v_y, const IntVect& iv) 
                   {
                       Real ry = (iv[1]+0.5)*h - 0.5;
                       Real rz = (iv[2]+0.5)*h - 0.5;
                       Real r = std::sqrt(ry*ry+rz*rz);
                       Real fac = std::exp(-(r*r/(0.16*0.16)));
                       v_y = v_y + 2.0*r*rz/r*fac;
                   });
               icomp = 2; // vz
               fab.ForEachIV(bx, icomp, ncomp, [=] (Real& v_z, const IntVect& iv) 
                   {
                       Real rz = (iv[2]+0.5)*h - 0.5;
                       Real ry = (iv[1]+0.5)*h - 0.5;
                       Real r = std::sqrt(ry*ry+rz*rz);
                       Real fac = std::exp(-(r*r/(0.16*0.16)));
                       v_z = v_z - 2.0*r*ry/r*fac;
                });
            } 
#endif
        }
    }
}
