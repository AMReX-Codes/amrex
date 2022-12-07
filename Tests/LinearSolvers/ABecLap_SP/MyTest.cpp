#include "MyTest.H"

#include <AMReX_ParmParse.H>

using namespace amrex;

MyTest::MyTest ()
{
    readParameters();
    initData();
}

void
MyTest::solve ()
{
    if (prob_type == 1) {
        if (single_precision) {
            solvePoisson<fMultiFab>();
        } else {
            solvePoisson<MultiFab>();
        }
    } else if (prob_type == 2) {
        if (single_precision) {
            solveABecLaplacian<fMultiFab>();
        } else {
            solveABecLaplacian<MultiFab>();
        }
    } else {
        amrex::Abort("Unknown prob_type");
    }
}

void
MyTest::readParameters ()
{
    ParmParse pp;
    pp.query("max_level", max_level);
    pp.query("ref_ratio", ref_ratio);
    pp.query("n_cell", n_cell);
    pp.query("max_grid_size", max_grid_size);

    pp.query("composite_solve", composite_solve);

    pp.query("prob_type", prob_type);

    pp.query("single_precision", single_precision);

    pp.query("verbose", verbose);
    pp.query("bottom_verbose", bottom_verbose);
    pp.query("max_iter", max_iter);
    pp.query("max_fmg_iter", max_fmg_iter);
    pp.query("linop_maxorder", linop_maxorder);
    pp.query("agglomeration", agglomeration);
    pp.query("consolidation", consolidation);
    pp.query("max_coarsening_level", max_coarsening_level);
}

void
MyTest::initData ()
{
    int nlevels = max_level + 1;
    geom.resize(nlevels);
    grids.resize(nlevels);
    dmap.resize(nlevels);

    solution.resize(nlevels);
    rhs.resize(nlevels);
    exact_solution.resize(nlevels);

    if (prob_type == 2 || prob_type == 3) {
        acoef.resize(nlevels);
        bcoef.resize(nlevels);
    }

    RealBox rb({AMREX_D_DECL(0.,0.,0.)}, {AMREX_D_DECL(1.,1.,1.)});
    Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0,0,0)};
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

    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        dmap[ilev].define(grids[ilev]);
        solution      [ilev].define(grids[ilev], dmap[ilev], 1, 1);
        rhs           [ilev].define(grids[ilev], dmap[ilev], 1, 0);
        exact_solution[ilev].define(grids[ilev], dmap[ilev], 1, 0);
        if (!acoef.empty()) {
            acoef[ilev].define(grids[ilev], dmap[ilev], 1, 0);
            bcoef[ilev].define(grids[ilev], dmap[ilev], 1, 1);
        }
    }

    if (prob_type == 1) {
        initProbPoisson();
    } else if (prob_type == 2) {
        initProbABecLaplacian();
    } else {
        amrex::Abort("Unknown prob_type "+std::to_string(prob_type));
    }
}
