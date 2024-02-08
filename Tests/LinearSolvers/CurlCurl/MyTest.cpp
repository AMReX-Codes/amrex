#include <AMReX_MLCurlCurl.H>

#include "MyTest.H"

#include <AMReX_MLMG.H>
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
    LPInfo info;
    info.setAgglomeration(agglomeration);
    info.setConsolidation(consolidation);
    info.setMaxCoarseningLevel(max_coarsening_level);

    MLCurlCurl mlcc({geom}, {grids}, {dmap}, info);

    mlcc.setDomainBC({AMREX_D_DECL(LinOpBCType::symmetry,
                                   LinOpBCType::Dirichlet,
                                   LinOpBCType::Periodic)},
                     {AMREX_D_DECL(LinOpBCType::Dirichlet,
                                   LinOpBCType::symmetry,
                                   LinOpBCType::Periodic)});

    mlcc.setScalars(alpha, beta);

    MLMGT<Array<MultiFab,3> > mlmg(mlcc);
    mlmg.setMaxIter(max_iter);
    mlmg.setVerbose(verbose);
    mlmg.setBottomVerbose(bottom_verbose);
    for (auto& mf : solution) {
        mf.setVal(Real(0));
    }
    mlmg.solve({&solution}, {&rhs}, Real(1.0e-10), Real(0));

    amrex::Print() << "  Number of cells: " << n_cell << std::endl;
    auto dvol = AMREX_D_TERM(geom.CellSize(0), *geom.CellSize(1), *geom.CellSize(2));
    Array<std::string,3> names{"Ex", "Ey", "Ez"};
    for (int idim = 0; idim < 3; ++idim) {
        MultiFab::Subtract(solution[idim], exact[idim], 0, 0, 1, 0);
        auto e0 = solution[idim].norminf();
        auto e1 = solution[idim].norm1(0,geom.periodicity());
        e1 *= dvol;
        auto e2 = solution[idim].norm2(0,geom.periodicity());
        e2 *= std::sqrt(dvol);
        amrex::Print() << "  " << names[idim] << " errors (max, L1, L2): "
                       << e0 << " " << e1 << " " << e2 << std::endl;
    }
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
    pp.query("agglomeration", agglomeration);
    pp.query("consolidation", consolidation);
    pp.query("max_coarsening_level", max_coarsening_level);

    pp.query("alpha_over_dx2", alpha_over_dx2);
    pp.query("beta", beta);
}

void
MyTest::initData ()
{
    RealBox rb({AMREX_D_DECL(0.,0.,0.)}, {AMREX_D_DECL(1.,1.,1.)});
    Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0,0,1)};
    Geometry::Setup(&rb, 0, is_periodic.data());
    Box domain(IntVect(0), IntVect(n_cell-1));
    geom.define(domain);

    const Real dx = geom.CellSize(0);
    alpha = alpha_over_dx2 * dx*dx;

    grids.define(domain);
    grids.maxSize(max_grid_size);
    dmap.define(grids);

    for (int idim = 0; idim < 3; ++idim) {
        IntVect itype(1);
#if (AMREX_SPACEDIM == 2)
        if (idim < AMREX_SPACEDIM)
#endif
        {
            itype[idim] = 0;
        }
        BoxArray const& ba = amrex::convert(grids, itype);
        solution[idim].define(ba,dmap,1,1);
        exact   [idim].define(ba,dmap,1,1);
        rhs     [idim].define(ba,dmap,1,0);
    }

    initProb();

    for (int idim = 0; idim < 3; ++idim) {
        exact[idim].LocalCopy(solution[idim], 0, 0, 1, IntVect(1));
    }
}
