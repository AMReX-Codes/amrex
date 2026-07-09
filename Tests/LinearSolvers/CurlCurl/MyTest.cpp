#include <AMReX_MLCurlCurl.H>

#include "MyTest.H"

#include <AMReX_GMRES.H>
#include <AMReX_GMRES_MLMG.H>
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
    Array<std::string,3> names{"Ex", "Ey", "Ez"};
    auto dvol = AMREX_D_TERM(geom.CellSize(0), *geom.CellSize(1), *geom.CellSize(2));

    LPInfo info;
    info.setAgglomeration(agglomeration);
    info.setConsolidation(consolidation);
    info.setMaxCoarseningLevel(max_coarsening_level);

    MLCurlCurl mlcc({geom}, {grids}, {dmap}, info, coord);

    mlcc.setDomainBC({AMREX_D_DECL(LinOpBCType::symmetry,
                                   LinOpBCType::Dirichlet,
                                   LinOpBCType::Periodic)},
                     {AMREX_D_DECL(LinOpBCType::Dirichlet,
                                   LinOpBCType::symmetry,
                                   LinOpBCType::Periodic)});

    mlcc.setScalars(alpha, beta);
    if (variable_beta) {
        Array<MultiFab,3> bcoef;
        for (int idim = 0; idim < 3; ++idim) {
            bcoef[idim].define(solution[idim].boxArray(),
                               solution[idim].DistributionMap(), 1, 0);
            bcoef[idim].setVal(beta);
        }
        mlcc.setBeta({Array<MultiFab const*,3>{bcoef.data(),
                                               bcoef.data()+1,
                                               bcoef.data()+2}});
    }
    if (variable_alpha) {
        mlcc.setAlpha({&alpha_node});
    }
    mlcc.prepareRHS({&rhs});

    if (use_pcg) { mlcc.setUsePCG(true); }

    using V = Array<MultiFab,3>;
    MLMGT<V> mlmg(mlcc);
    mlmg.setMaxIter(max_iter);
    mlmg.setVerbose(verbose);
    mlmg.setBottomVerbose(bottom_verbose);

    {
        Array<MultiFab,3> ax_exact;
        Array<MultiFab,3> residual;
        for (int idim = 0; idim < 3; ++idim) {
            ax_exact[idim].define(rhs[idim].boxArray(), rhs[idim].DistributionMap(), 1, 0);
            residual[idim].define(rhs[idim].boxArray(), rhs[idim].DistributionMap(), 1, 0);
            MultiFab::Copy(residual[idim], rhs[idim], 0, 0, 1, 0);
        }
        mlmg.apply({&ax_exact}, {&exact});
        amrex::Print() << "  Residual norms (b - A x_exact), expect O(dx^2):\n";
        for (int idim = 0; idim < 3; ++idim) {
            MultiFab::Subtract(residual[idim], ax_exact[idim], 0, 0, 1, 0);
        }
        mlcc.setDirichletNodesToZero(0,0,residual);
        for (int idim = 0; idim < 3; ++idim) {
            auto r0 = residual[idim].norminf();
            auto r1 = residual[idim].norm1(0, geom.periodicity()) * dvol;
            auto r2 = residual[idim].norm2(0, geom.periodicity()) * std::sqrt(dvol);
            amrex::Print() << "    " << names[idim] << ": "
                           << r0 << " " << r1 << " " << r2 << '\n';
        }
    }

    for (auto& mf : solution) {
        mf.setVal(Real(0));
    }

    Real tol_rel;
    if constexpr (std::is_same_v<double,Real>) {
        tol_rel = Real(1.0e-10);
    } else {
        tol_rel = Real(1.0e-4);
    }
    auto tol_abs = Real(0.0);

    if (use_gmres)
    {
        GMRESMLMGT<V> gmsolver(mlmg);
        gmsolver.usePrecond(gmres_use_precond);
        gmsolver.setPrecondNumIters(gmres_precond_niters);
        gmsolver.setVerbose(verbose);
        gmsolver.solve(solution, rhs, tol_rel, tol_abs);
    }
    else
    {
        mlmg.solve({&solution}, {&rhs}, tol_rel, tol_abs);
    }

    amrex::Print() << "  Number of cells: " << n_cell << '\n';
    for (int idim = 0; idim < 3; ++idim) {
        MultiFab::Subtract(solution[idim], exact[idim], 0, 0, 1, 0);
        auto e0 = solution[idim].norminf();
        auto e1 = solution[idim].norm1(0,geom.periodicity());
        e1 *= dvol;
        auto e2 = solution[idim].norm2(0,geom.periodicity());
        e2 *= std::sqrt(dvol);
        amrex::Print() << "  " << names[idim] << " errors (max, L1, L2): "
                       << e0 << " " << e1 << " " << e2 << '\n';
    }
}

void
MyTest::readParameters ()
{
    ParmParse pp;
    pp.query("coord", coord);
#if (AMREX_SPACEDIM == 3)
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(coord == 0, "Only Cartesian is supported in 3D");
#elif (AMREX_SPACEDIM == 2)
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(coord != 2, "Spherical coordinates not supported in 2D");
#endif

    pp.query("n_cell", n_cell);
    pp.query("max_grid_size", max_grid_size);

    pp.query("verbose", verbose);
    pp.query("bottom_verbose", bottom_verbose);
    pp.query("max_iter", max_iter);
    pp.query("agglomeration", agglomeration);
    pp.query("consolidation", consolidation);
    pp.query("max_coarsening_level", max_coarsening_level);

    pp.query("use_pcg", use_pcg);
    pp.query("use_gmres", use_gmres);
    pp.query("gmres_use_precond", gmres_use_precond);
    pp.query("gmres_precond_niters", gmres_precond_niters);

    pp.query("beta_factor", beta_factor);
    pp.query("alpha", alpha);
    pp.query("variable_beta", variable_beta);
    pp.query("variable_alpha", variable_alpha);
}

void
MyTest::initData ()
{
    RealBox rb(AMREX_D_DECL(0.,0.,0.), AMREX_D_DECL(1.,1.,1.));
    Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0,0,1)};
    Geometry::Setup(&rb, 0, is_periodic.data());
    Box domain(IntVect(0), IntVect(n_cell-1));
    geom.define(domain);

    const Real dx = geom.CellSize(0);
    beta = beta_factor * alpha/(dx*dx);

    grids.define(domain);
    grids.maxSize(max_grid_size);
    dmap.define(grids);

    for (int idim = 0; idim < 3; ++idim) {
        IntVect itype(1);
#if (AMREX_SPACEDIM < 3)
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
    if (variable_alpha) {
        alpha_node.define(amrex::convert(grids,IntVect(1)), dmap,1,0);
    }

    initProb();

    for (int idim = 0; idim < 3; ++idim) {
        exact[idim].LocalCopy(solution[idim], 0, 0, 1, IntVect(1));
    }
}
