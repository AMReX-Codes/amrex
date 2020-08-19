#include "MyTest.H"

#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>

using namespace amrex;

MyTest::MyTest ()
{
    readParameters();

    initGrids();

    initData();
}

//
// Solve L(phi) = rhs
//
void
MyTest::solve ()
{
    std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_lobc;
    std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_hibc;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        mlmg_lobc[idim] = LinOpBCType::Dirichlet;
        mlmg_hibc[idim] = LinOpBCType::Dirichlet;
    }

    LPInfo info;
    info.setMaxCoarseningLevel(max_coarsening_level);

    MLNodeLaplacian mlndlap({geom}, {grids}, {dmap}, info);

    mlndlap.setDomainBC(mlmg_lobc, mlmg_hibc);

    mlndlap.setOversetMask(0, dmask);

    {
        MultiFab sigma(grids, dmap, 1, 0);
        sigma.setVal(1.0);
        mlndlap.setSigma(0, sigma);
    }

    MLMG mlmg(mlndlap);
    mlmg.setVerbose(verbose);
    mlmg.setBottomVerbose(bottom_verbose);

    const auto domain = geom.Domain();
    for (MFIter mfi(phi); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();
        Array4<Real> const& pa = phi.array(mfi);
        Array4<Real> const& ra = rhs.array(mfi);
        Array4<int const> const& ma = dmask.array(mfi);
        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (domain.strictly_contains(IntVect(AMREX_D_DECL(i,j,k))) and
                ma(i,j,k) == 1)
            {  // Let's set phi = 0 for unknown nodes
                pa(i,j,k) = 0.0;
            }
            if (ma(i,j,k) == 0) {  // Let's set rhs = 0 for masked out known nodes
                ra(i,j,k) = 0.0;
            }
        });
    }

    Real mlmg_err = mlmg.solve({&phi}, {&rhs}, 1.e-11, 0.0);

#if 1
    VisMF::Write(phi, "phi");
    VisMF::Write(rhs, "rhs");
#endif
}

void
MyTest::writePlotfile ()
{
//    amrex::Print() << "MyTest::writePlotfile todo" << std::endl;
}

void
MyTest::readParameters ()
{
    ParmParse pp;
    pp.query("n_cell", n_cell);
    pp.query("max_grid_size", max_grid_size);

    pp.query("plot_file", plot_file_name);

    pp.query("verbose", verbose);
    pp.query("bottom_verbose", bottom_verbose);
    pp.query("max_coarsening_level", max_coarsening_level);
}

void
MyTest::initGrids ()
{
    RealBox rb({AMREX_D_DECL(0.,0.,0.)}, {AMREX_D_DECL(1.,1.,1.)});
    std::array<int,AMREX_SPACEDIM> isperiodic{AMREX_D_DECL(0,0,0)};
    Geometry::Setup(&rb, 0, isperiodic.data());
    Box domain(IntVect{AMREX_D_DECL(0,0,0)}, IntVect{AMREX_D_DECL(n_cell-1,n_cell-1,n_cell-1)});
    geom.define(domain, rb, CoordSys::cartesian, isperiodic);

    grids.define(domain);
    grids.maxSize(max_grid_size);
}

void
MyTest::initData ()
{
    dmap.define(grids);
    const BoxArray& nba = amrex::convert(grids,IntVect::TheNodeVector());
    phi.define(nba, dmap, 1, 0);
    rhs.define(nba, dmap, 1, 0);
    dmask.define(nba, dmap, 1, 0);

    const auto dx = geom.CellSizeArray();
    static_assert(AMREX_SPACEDIM == 2, "3D todo");
    for (MFIter mfi(phi); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();
        Array4<Real> const& pa = phi.array(mfi);
        Array4<Real> const& ra = rhs.array(mfi);
        Array4<int> const& ma = dmask.array(mfi);
        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real x = i*dx[0] - 0.5;
            Real y = j*dx[1] - 0.5;
            Real theta = std::atan2(x,y) + 0.5*3.1415926535897932;
            Real r2 = x*x + y*y;
            pa(i,j,k) =  r2*r2*std::cos(3.0*theta);
            ra(i,j,k) = 7.0*r2*std::cos(3.0*theta);
            Real r = std::sqrt(r2);
            if (r < (0.3 + 0.15*std::cos(6.*theta))) {
                ma(i,j,k) = 0; // masked out known nodes
            } else {
                ma(i,j,k) = 1;
            }
        });
    }
}
