#include "MyTest.H"
#include "MyTest_F.H"

#include <AMReX_MLNodeLaplacian.H>
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
    MLNodeLaplacian linop({geom}, {grids}, {dmap});

    linop.setDomainBC({AMREX_D_DECL(LinOpBCType::Periodic,
                                    LinOpBCType::Periodic,
                                    LinOpBCType::Periodic)},
                      {AMREX_D_DECL(LinOpBCType::Periodic,
                                    LinOpBCType::Periodic,
                                    LinOpBCType::Periodic)});

    linop.setLevelBC(0, nullptr);

    linop.setSigma(0, sigma);

    MLMG mlmg(linop);
    mlmg.setMaxIter(max_iter);
    mlmg.setMaxFmgIter(max_fmg_iter);
    mlmg.setVerbose(verbose);
    mlmg.setBottomVerbose(bottom_verbose);

    mlmg.solve({&solution}, {&rhs}, reltol, 0.0);
}

void
MyTest::compute_norms () const
{
    MultiFab error(solution.boxArray(), solution.DistributionMap(), 1, 0);
    MultiFab::Copy(error, solution, 0, 0, 1, 0);
    MultiFab::Subtract(error, exact_solution, 0, 0, 1, 0);

    auto mask = error.OwnerMask(geom.periodicity());

    // compute average
    Real avg;
    {
        MultiFab one(error.boxArray(), error.DistributionMap(), 1, 0);
        one.setVal(1.0);
        Real r1 = MultiFab::Dot(*mask,one,0,one,0,1,0);
        Real r2 = MultiFab::Dot(*mask,one,0,error,0,1,0);
        avg = r2 / r1;
    }

    error.plus(-avg,0,1,0); // so that the sum is zero

    amrex::Print() << "max-norm: " << error.norm0(*mask, 0, 0) << "\n";
    const Real* dx = geom.CellSize();
    Real dvol = AMREX_D_TERM(dx[0], *dx[1], *dx[2]);
    amrex::Print() << "1-norm  : " << error.norm1(0, geom.periodicity())*dvol << "\n";
}

void
MyTest::readParameters ()
{
    ParmParse pp;
    pp.query("n_cell", n_cell);
    pp.query("max_grid_size", max_grid_size);

    pp.query("verbose", verbose);
    pp.query("bottom_verbose", bottom_verbose);
    pp.query("max_iter", max_iter);
    pp.query("max_fmg_iter", max_fmg_iter);
    pp.query("reltol", reltol);
}

void
MyTest::initData ()
{
    RealBox rb({AMREX_D_DECL(0.,0.,0.)}, {AMREX_D_DECL(1.,1.,1.)});
    Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1,1,1)};
    Geometry::Setup(&rb, 0, is_periodic.data());
    Box domain0(IntVect{AMREX_D_DECL(0,0,0)}, IntVect{AMREX_D_DECL(n_cell-1,n_cell-1,n_cell-1)});
    geom.define(domain0);

    grids.define(domain0);
    grids.maxSize(max_grid_size);

    dmap.define(grids);

    solution.define      (amrex::convert(grids,IntVect::TheUnitVector()), dmap, 1, 0);
    rhs.define           (amrex::convert(grids,IntVect::TheUnitVector()), dmap, 1, 0);
    exact_solution.define(amrex::convert(grids,IntVect::TheUnitVector()), dmap, 1, 0);
    sigma.define         (               grids,                           dmap, 1, 0);

    solution.setVal(0.0);
    sigma.setVal(1.0);

    const Real* dx = geom.CellSize();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(rhs,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        init_prob(BL_TO_FORTRAN_BOX(bx),
                  BL_TO_FORTRAN_ANYD(rhs[mfi]),
                  BL_TO_FORTRAN_ANYD(exact_solution[mfi]),
                  dx);
    }
}

