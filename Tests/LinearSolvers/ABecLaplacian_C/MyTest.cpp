#include "MyTest.H"

#include <AMReX_GMRES.H>
#include <AMReX_GMRES_MLMG.H>
#include <AMReX_MLNodeABecLaplacian.H>
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLPoisson.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>

#ifdef AMREX_USE_HYPRE
#include <AMReX_HypreMLABecLap.H>
#endif

using namespace amrex;

MyTest::MyTest ()
{
    readParameters();
    initData();
}

void
MyTest::solve ()
{
#ifdef AMREX_USE_HYPRE
    if (use_mlhypre) {
        solveMLHypre();
        return;
    }
#endif

    if (prob_type == 1) {
        solvePoisson();
    } else if (prob_type == 2) {
        if (use_gmres) {
            solveABecLaplacianGMRES();
        } else {
            solveABecLaplacian();
        }
    } else if (prob_type == 3) {
        solveABecLaplacianInhomNeumann();
    } else if (prob_type == 4) {
        solveNodeABecLaplacian();
    } else {
        amrex::Abort("Unknown prob_type");
    }
}

void
MyTest::solvePoisson ()
{
    LPInfo info;
    info.setAgglomeration(agglomeration);
    info.setConsolidation(consolidation);
    info.setMaxCoarseningLevel(max_coarsening_level);

    const auto tol_rel = Real(1.e-10);
    const auto tol_abs = Real(0.0);

    const auto nlevels = static_cast<int>(geom.size());

    if (composite_solve)
    {
        MLPoisson mlpoisson(geom, grids, dmap, info);

        mlpoisson.setGaussSeidel(use_gauss_seidel);

        mlpoisson.setMaxOrder(linop_maxorder);

        // This is a 3d problem with Dirichlet BC
        mlpoisson.setDomainBC({AMREX_D_DECL(LinOpBCType::Dirichlet,
                                            LinOpBCType::Dirichlet,
                                            LinOpBCType::Dirichlet)},
                              {AMREX_D_DECL(LinOpBCType::Dirichlet,
                                            LinOpBCType::Dirichlet,
                                            LinOpBCType::Dirichlet)});

        for (int ilev = 0; ilev < nlevels; ++ilev)
        {
            mlpoisson.setLevelBC(ilev, &solution[ilev]);
        }

        MLMG mlmg(mlpoisson);
        mlmg.setMaxIter(max_iter);
        mlmg.setMaxFmgIter(max_fmg_iter);
        mlmg.setVerbose(verbose);
        mlmg.setBottomVerbose(bottom_verbose);
#ifdef AMREX_USE_HYPRE
        if (use_hypre) {
            mlmg.setBottomSolver(MLMG::BottomSolver::hypre);
            mlmg.setHypreInterface(hypre_interface);
        }
#endif
#ifdef AMREX_USE_PETSC
        if (use_petsc) {
            mlmg.setBottomSolver(MLMG::BottomSolver::petsc);
        }
#endif

        mlmg.solve(GetVecOfPtrs(solution), GetVecOfConstPtrs(rhs), tol_rel, tol_abs);
    }
    else
    {
        for (int ilev = 0; ilev < nlevels; ++ilev)
        {
            MLPoisson mlpoisson({geom[ilev]}, {grids[ilev]}, {dmap[ilev]}, info);

            mlpoisson.setGaussSeidel(use_gauss_seidel);

            mlpoisson.setMaxOrder(linop_maxorder);

            // This is a 3d problem with Dirichlet BC
            mlpoisson.setDomainBC({AMREX_D_DECL(LinOpBCType::Dirichlet,
                                                LinOpBCType::Dirichlet,
                                                LinOpBCType::Dirichlet)},
                                  {AMREX_D_DECL(LinOpBCType::Dirichlet,
                                                LinOpBCType::Dirichlet,
                                                LinOpBCType::Dirichlet)});

            if (ilev > 0) {
                mlpoisson.setCoarseFineBC(&solution[ilev-1], ref_ratio);
            }

            mlpoisson.setLevelBC(0, &solution[ilev]);

            MLMG mlmg(mlpoisson);
            mlmg.setMaxIter(max_iter);
            mlmg.setMaxFmgIter(max_fmg_iter);
            mlmg.setVerbose(verbose);
            mlmg.setBottomVerbose(bottom_verbose);
#ifdef AMREX_USE_HYPRE
            if (use_hypre) {
                mlmg.setBottomSolver(MLMG::BottomSolver::hypre);
                mlmg.setHypreInterface(hypre_interface);
            }
#endif
#ifdef AMREX_USE_PETSC
            if (use_petsc) {
                mlmg.setBottomSolver(MLMG::BottomSolver::petsc);
            }
#endif

            mlmg.solve({&solution[ilev]}, {&rhs[ilev]}, tol_rel, tol_abs);
        }
    }
}

void
MyTest::solveABecLaplacian ()
{
    LPInfo info;
    info.setAgglomeration(agglomeration);
    info.setConsolidation(consolidation);
    info.setSemicoarsening(semicoarsening);
    info.setMaxCoarseningLevel(max_coarsening_level);
    info.setMaxSemicoarseningLevel(max_semicoarsening_level);

    const auto tol_rel = Real(1.e-10);
    const auto tol_abs = Real(0.0);

    const auto nlevels = static_cast<int>(geom.size());

    if (composite_solve)
    {

        MLABecLaplacian mlabec(geom, grids, dmap, info);

        mlabec.setGaussSeidel(use_gauss_seidel);

        mlabec.setMaxOrder(linop_maxorder);

        mlabec.setDomainBC({AMREX_D_DECL(LinOpBCType::Dirichlet,
                                         LinOpBCType::Neumann,
                                         LinOpBCType::Neumann)},
                           {AMREX_D_DECL(LinOpBCType::Neumann,
                                         LinOpBCType::Dirichlet,
                                         LinOpBCType::Neumann)});

        for (int ilev = 0; ilev < nlevels; ++ilev)
        {
            mlabec.setLevelBC(ilev, &solution[ilev]);
        }

        mlabec.setScalars(ascalar, bscalar);

        for (int ilev = 0; ilev < nlevels; ++ilev)
        {
            mlabec.setACoeffs(ilev, acoef[ilev]);

            Array<MultiFab,AMREX_SPACEDIM> face_bcoef;
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
            {
                const BoxArray& ba = amrex::convert(bcoef[ilev].boxArray(),
                                                    IntVect::TheDimensionVector(idim));
                face_bcoef[idim].define(ba, bcoef[ilev].DistributionMap(), 1, 0);
            }
            amrex::average_cellcenter_to_face(GetArrOfPtrs(face_bcoef),
                                              bcoef[ilev], geom[ilev]);
            mlabec.setBCoeffs(ilev, amrex::GetArrOfConstPtrs(face_bcoef));
        }

        MLMG mlmg(mlabec);
        mlmg.setMaxIter(max_iter);
        mlmg.setMaxFmgIter(max_fmg_iter);
        mlmg.setVerbose(verbose);
        mlmg.setBottomVerbose(bottom_verbose);
#ifdef AMREX_USE_HYPRE
        if (use_hypre) {
            mlmg.setBottomSolver(MLMG::BottomSolver::hypre);
            mlmg.setHypreInterface(hypre_interface);
        }
#endif
#ifdef AMREX_USE_PETSC
        if (use_petsc) {
            mlmg.setBottomSolver(MLMG::BottomSolver::petsc);
        }
#endif

        mlmg.solve(GetVecOfPtrs(solution), GetVecOfConstPtrs(rhs), tol_rel, tol_abs);
    }
    else
    {
        for (int ilev = 0; ilev < nlevels; ++ilev)
        {
            MLABecLaplacian mlabec({geom[ilev]}, {grids[ilev]}, {dmap[ilev]}, info);

            mlabec.setGaussSeidel(use_gauss_seidel);

            mlabec.setMaxOrder(linop_maxorder);

            mlabec.setDomainBC({AMREX_D_DECL(LinOpBCType::Dirichlet,
                                             LinOpBCType::Neumann,
                                             LinOpBCType::Neumann)},
                               {AMREX_D_DECL(LinOpBCType::Neumann,
                                             LinOpBCType::Dirichlet,
                                             LinOpBCType::Neumann)});

            if (ilev > 0) {
                mlabec.setCoarseFineBC(&solution[ilev-1], ref_ratio);
            }

            mlabec.setLevelBC(0, &solution[ilev]);

            mlabec.setScalars(ascalar, bscalar);

            mlabec.setACoeffs(0, acoef[ilev]);

            Array<MultiFab,AMREX_SPACEDIM> face_bcoef;
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
            {
                const BoxArray& ba = amrex::convert(bcoef[ilev].boxArray(),
                                                    IntVect::TheDimensionVector(idim));
                face_bcoef[idim].define(ba, bcoef[ilev].DistributionMap(), 1, 0);
            }
            amrex::average_cellcenter_to_face(GetArrOfPtrs(face_bcoef),
                                              bcoef[ilev], geom[ilev]);
            mlabec.setBCoeffs(0, amrex::GetArrOfConstPtrs(face_bcoef));

            MLMG mlmg(mlabec);
            mlmg.setMaxIter(max_iter);
            mlmg.setMaxFmgIter(max_fmg_iter);
            mlmg.setVerbose(verbose);
            mlmg.setBottomVerbose(bottom_verbose);
#ifdef AMREX_USE_HYPRE
            if (use_hypre) {
                mlmg.setBottomSolver(MLMG::BottomSolver::hypre);
                mlmg.setHypreInterface(hypre_interface);
            }
#endif
#ifdef AMREX_USE_PETSC
            if (use_petsc) {
                mlmg.setBottomSolver(MLMG::BottomSolver::petsc);
            }
#endif

            mlmg.solve({&solution[ilev]}, {&rhs[ilev]}, tol_rel, tol_abs);
        }
    }
}

void
MyTest::solveABecLaplacianInhomNeumann ()
{
    LPInfo info;
    info.setAgglomeration(agglomeration);
    info.setConsolidation(consolidation);
    info.setMaxCoarseningLevel(max_coarsening_level);

    const auto tol_rel = Real(1.e-10);
    const auto tol_abs = Real(0.0);

    const auto nlevels = static_cast<int>(geom.size());

    if (composite_solve)
    {

        MLABecLaplacian mlabec(geom, grids, dmap, info);

        mlabec.setGaussSeidel(use_gauss_seidel);

        mlabec.setMaxOrder(linop_maxorder);

        // This is a 3d problem with inhomogeneous Neumann BC
        mlabec.setDomainBC({AMREX_D_DECL(LinOpBCType::inhomogNeumann,
                                         LinOpBCType::inhomogNeumann,
                                         LinOpBCType::inhomogNeumann)},
                           {AMREX_D_DECL(LinOpBCType::inhomogNeumann,
                                         LinOpBCType::inhomogNeumann,
                                         LinOpBCType::inhomogNeumann)});

        for (int ilev = 0; ilev < nlevels; ++ilev)
        {
            // for problem with inhomogeneous Neumann BC, we need to pass boundary values
            mlabec.setLevelBC(ilev, &(solution[ilev]));
        }

        mlabec.setScalars(ascalar, bscalar);

        for (int ilev = 0; ilev < nlevels; ++ilev)
        {
            mlabec.setACoeffs(ilev, acoef[ilev]);

            Array<MultiFab,AMREX_SPACEDIM> face_bcoef;
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
            {
                const BoxArray& ba = amrex::convert(bcoef[ilev].boxArray(),
                                                    IntVect::TheDimensionVector(idim));
                face_bcoef[idim].define(ba, bcoef[ilev].DistributionMap(), 1, 0);
            }
            amrex::average_cellcenter_to_face(GetArrOfPtrs(face_bcoef),
                                              bcoef[ilev], geom[ilev]);
            mlabec.setBCoeffs(ilev, amrex::GetArrOfConstPtrs(face_bcoef));
        }

        MLMG mlmg(mlabec);
        mlmg.setMaxIter(max_iter);
        mlmg.setMaxFmgIter(max_fmg_iter);
        mlmg.setVerbose(verbose);
        mlmg.setBottomVerbose(bottom_verbose);
#ifdef AMREX_USE_HYPRE
        if (use_hypre) {
            mlmg.setBottomSolver(MLMG::BottomSolver::hypre);
            mlmg.setHypreInterface(hypre_interface);
        }
#endif
#ifdef AMREX_USE_PETSC
        if (use_petsc) {
            mlmg.setBottomSolver(MLMG::BottomSolver::petsc);
        }
#endif

        mlmg.solve(GetVecOfPtrs(solution), GetVecOfConstPtrs(rhs), tol_rel, tol_abs);
    }
    else
    {
        for (int ilev = 0; ilev < nlevels; ++ilev)
        {
            MLABecLaplacian mlabec({geom[ilev]}, {grids[ilev]}, {dmap[ilev]}, info);

            mlabec.setGaussSeidel(use_gauss_seidel);

            mlabec.setMaxOrder(linop_maxorder);

            // This is a 3d problem with inhomogeneous Neumann BC
            mlabec.setDomainBC({AMREX_D_DECL(LinOpBCType::inhomogNeumann,
                                             LinOpBCType::inhomogNeumann,
                                             LinOpBCType::inhomogNeumann)},
                               {AMREX_D_DECL(LinOpBCType::inhomogNeumann,
                                             LinOpBCType::inhomogNeumann,
                                             LinOpBCType::inhomogNeumann)});

            if (ilev > 0) {
                mlabec.setCoarseFineBC(&solution[ilev-1], ref_ratio);
            }

            // for problem with inhomogeneous Neumann BC, we need to pass boundary values
            mlabec.setLevelBC(0, &solution[ilev]);

            mlabec.setScalars(ascalar, bscalar);

            mlabec.setACoeffs(0, acoef[ilev]);

            Array<MultiFab,AMREX_SPACEDIM> face_bcoef;
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
            {
                const BoxArray& ba = amrex::convert(bcoef[ilev].boxArray(),
                                                    IntVect::TheDimensionVector(idim));
                face_bcoef[idim].define(ba, bcoef[ilev].DistributionMap(), 1, 0);
            }
            amrex::average_cellcenter_to_face(GetArrOfPtrs(face_bcoef),
                                              bcoef[ilev], geom[ilev]);
            mlabec.setBCoeffs(0, amrex::GetArrOfConstPtrs(face_bcoef));

            MLMG mlmg(mlabec);
            mlmg.setMaxIter(max_iter);
            mlmg.setMaxFmgIter(max_fmg_iter);
            mlmg.setVerbose(verbose);
            mlmg.setBottomVerbose(bottom_verbose);
#ifdef AMREX_USE_HYPRE
            if (use_hypre) {
                mlmg.setBottomSolver(MLMG::BottomSolver::hypre);
                mlmg.setHypreInterface(hypre_interface);
            }
#endif
#ifdef AMREX_USE_PETSC
            if (use_petsc) {
                mlmg.setBottomSolver(MLMG::BottomSolver::petsc);
            }
#endif

            mlmg.solve({&solution[ilev]}, {&rhs[ilev]}, tol_rel, tol_abs);
        }
    }
}

void
MyTest::solveNodeABecLaplacian ()
{
    LPInfo info;
    info.setAgglomeration(agglomeration);
    info.setConsolidation(consolidation);
    info.setMaxCoarseningLevel(max_coarsening_level);

    const auto tol_rel = Real(1.e-10);
    const auto tol_abs = Real(0.0);

    const auto nlevels = static_cast<int>(geom.size());

    if (composite_solve && nlevels > 1)
    {
        amrex::Abort("solveNodeABecLaplacian: TODO composite_solve");
    }
    else
    {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(nlevels == 1, "solveNodeABecLaplacian: nlevels > 1 TODO");
        for (int ilev = 0; ilev < nlevels; ++ilev)
        {
            MLNodeABecLaplacian mlndabec({geom[ilev]}, {grids[ilev]}, {dmap[ilev]},
                                         info);

            mlndabec.setDomainBC({AMREX_D_DECL(LinOpBCType::Dirichlet,
                                               LinOpBCType::Neumann,
                                               LinOpBCType::Dirichlet)},
                                 {AMREX_D_DECL(LinOpBCType::Neumann,
                                               LinOpBCType::Dirichlet,
                                               LinOpBCType::Dirichlet)});

            mlndabec.setScalars(ascalar, bscalar);

            mlndabec.setACoeffs(0, acoef[ilev]);
            mlndabec.setBCoeffs(0, bcoef[ilev]);

            MLMG mlmg(mlndabec);
            mlmg.setMaxIter(max_iter);
            mlmg.setMaxFmgIter(max_fmg_iter);
            mlmg.setVerbose(verbose);
            mlmg.setBottomVerbose(bottom_verbose);

            mlmg.solve({&solution[ilev]}, {&rhs[ilev]}, tol_rel, tol_abs);
        }
    }
}

void
MyTest::solveABecLaplacianGMRES ()
{
    LPInfo info;
    info.setAgglomeration(agglomeration);
    info.setConsolidation(consolidation);
    info.setSemicoarsening(semicoarsening);
    info.setMaxCoarseningLevel(max_coarsening_level);
    info.setMaxSemicoarseningLevel(max_semicoarsening_level);

    const auto tol_rel = Real(1.e-10);
    const auto tol_abs = Real(0.0);

    const auto nlevels = static_cast<int>(geom.size());

    if (composite_solve)
    {
        MLABecLaplacian mlabec(geom, grids, dmap, info);

        mlabec.setMaxOrder(linop_maxorder);

        mlabec.setDomainBC({AMREX_D_DECL(LinOpBCType::Dirichlet,
                                         LinOpBCType::Neumann,
                                         LinOpBCType::Neumann)},
                           {AMREX_D_DECL(LinOpBCType::Neumann,
                                         LinOpBCType::Dirichlet,
                                         LinOpBCType::Neumann)});

        for (int ilev = 0; ilev < nlevels; ++ilev)
        {
            mlabec.setLevelBC(ilev, &solution[ilev]);
        }

        mlabec.setScalars(ascalar, bscalar);

        for (int ilev = 0; ilev < nlevels; ++ilev)
        {
            mlabec.setACoeffs(ilev, acoef[ilev]);

            Array<MultiFab,AMREX_SPACEDIM> face_bcoef;
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
            {
                const BoxArray& ba = amrex::convert(bcoef[ilev].boxArray(),
                                                    IntVect::TheDimensionVector(idim));
                face_bcoef[idim].define(ba, bcoef[ilev].DistributionMap(), 1, 0);
            }
            amrex::average_cellcenter_to_face(GetArrOfPtrs(face_bcoef),
                                              bcoef[ilev], geom[ilev]);
            mlabec.setBCoeffs(ilev, amrex::GetArrOfConstPtrs(face_bcoef));
        }

        MLMG mlmg(mlabec);
        GMRESMLMGT<MultiFab> gmsolver(mlmg);
        gmsolver.usePrecond(true);
        gmsolver.setVerbose(verbose);
        gmsolver.solve(GetVecOfPtrs(solution), GetVecOfConstPtrs(rhs), tol_rel, tol_abs);

        if (verbose) {
            Vector<MultiFab> res(nlevels);
            for (int ilev = 0; ilev < nlevels; ++ilev) {
                res[ilev].define(rhs[ilev].boxArray(), rhs[ilev].DistributionMap(), 1, 0);
            }
            mlmg.apply(GetVecOfPtrs(res), GetVecOfPtrs(solution)); // res = L(sol)
            for (int ilev = 0; ilev < nlevels; ++ilev) {
                MultiFab::Subtract(res[ilev], rhs[ilev], 0, 0, 1, 0);
                amrex::Print() << "Final residual at level " << ilev << " = "
                               << res[ilev].norminf(0) << " " << res[ilev].norm1(0) << " "
                               << res[ilev].norm2(0) << '\n';
            }
        }
    }
    else
    {
        for (int ilev = 0; ilev < nlevels; ++ilev)
        {
            MLABecLaplacian mlabec({geom[ilev]}, {grids[ilev]}, {dmap[ilev]}, info);

            mlabec.setMaxOrder(linop_maxorder);

            mlabec.setDomainBC({AMREX_D_DECL(LinOpBCType::Dirichlet,
                                             LinOpBCType::Neumann,
                                             LinOpBCType::Neumann)},
                               {AMREX_D_DECL(LinOpBCType::Neumann,
                                             LinOpBCType::Dirichlet,
                                             LinOpBCType::Neumann)});

            if (ilev > 0) {
                mlabec.setCoarseFineBC(&solution[ilev-1], ref_ratio);
            }

            // for problem with pure homogeneous Neumann BC, we could pass a nullptr
            mlabec.setLevelBC(0, &solution[ilev]);

            mlabec.setScalars(ascalar, bscalar);

            mlabec.setACoeffs(0, acoef[ilev]);

            Array<MultiFab,AMREX_SPACEDIM> face_bcoef;
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
            {
                const BoxArray& ba = amrex::convert(bcoef[ilev].boxArray(),
                                                    IntVect::TheDimensionVector(idim));
                face_bcoef[idim].define(ba, bcoef[ilev].DistributionMap(), 1, 0);
            }
            amrex::average_cellcenter_to_face(GetArrOfPtrs(face_bcoef),
                                              bcoef[ilev], geom[ilev]);
            mlabec.setBCoeffs(0, amrex::GetArrOfConstPtrs(face_bcoef));

            MLMG mlmg(mlabec);
            GMRESMLMGT gmsolver(mlmg);
            gmsolver.usePrecond(true);
            gmsolver.setVerbose(verbose);
            gmsolver.solve(solution[ilev], rhs[ilev], tol_rel, tol_abs);

            if (verbose) {
                MultiFab res(rhs[ilev].boxArray(), rhs[ilev].DistributionMap(), 1, 0);
                mlmg.apply({&res}, {&solution[ilev]}); // res = L(sol)
                MultiFab::Subtract(res, rhs[ilev], 0, 0, 1, 0); // now res = L(sol) - rhs
                amrex::Print() << "Final residual = " << res.norminf(0)
                               << " " << res.norm1(0) << " " << res.norm2(0) << '\n';
            }
        }
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

    pp.query("use_mlhypre", use_mlhypre);
    pp.query("use_hypre_ssamg", use_hypre_ssamg);

    pp.query("prob_type", prob_type);

    pp.query("verbose", verbose);
    pp.query("bottom_verbose", bottom_verbose);
    pp.query("max_iter", max_iter);
    pp.query("max_fmg_iter", max_fmg_iter);
    pp.query("linop_maxorder", linop_maxorder);
    pp.query("agglomeration", agglomeration);
    pp.query("consolidation", consolidation);
    pp.query("semicoarsening", semicoarsening);
    pp.query("max_coarsening_level", max_coarsening_level);
    pp.query("max_semicoarsening_level", max_semicoarsening_level);

    pp.query("use_gauss_seidel", use_gauss_seidel);

    pp.query("use_gmres", use_gmres);
    AMREX_ALWAYS_ASSERT(use_gmres == false || prob_type == 2);

#ifdef AMREX_USE_HYPRE
    pp.query("use_hypre", use_hypre);
    pp.query("hypre_interface", hypre_interface_i);
    if (hypre_interface_i == 1) {
        hypre_interface = Hypre::Interface::structed;
    } else if (hypre_interface_i == 2) {
        hypre_interface = Hypre::Interface::semi_structed;
    } else {
        hypre_interface = Hypre::Interface::ij;
    }
#endif
#ifdef AMREX_USE_PETSC
    pp.query("use_petsc", use_petsc);
#endif
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!(use_hypre && use_petsc),
                                     "use_hypre & use_petsc cannot be both true");
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

    if (prob_type == 2 || prob_type == 3 || prob_type == 4) {
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
        BoxArray ba = grids[ilev];
        if (prob_type == 4) {
            ba.surroundingNodes();
        }
        solution      [ilev].define(ba, dmap[ilev], 1, 1);
        rhs           [ilev].define(ba, dmap[ilev], 1, 0);
        exact_solution[ilev].define(ba, dmap[ilev], 1, 0);
        if (!acoef.empty()) {
            acoef[ilev].define(ba         , dmap[ilev], 1, 0);
            const int ngb = (prob_type == 4) ? 0 : 1;
            bcoef[ilev].define(grids[ilev], dmap[ilev], 1, ngb);
        }
    }

    if (prob_type == 1) {
        initProbPoisson();
    } else if (prob_type == 2) {
        initProbABecLaplacian();
    } else if (prob_type == 3) {
        initProbABecLaplacianInhomNeumann();
    } else if (prob_type == 4) {
        initProbNodeABecLaplacian();
    } else {
        amrex::Abort("Unknown prob_type "+std::to_string(prob_type));
    }
}

#ifdef AMREX_USE_HYPRE
void
MyTest::solveMLHypre ()
{
    const auto tol_rel = Real(1.e-10);
    const auto tol_abs = Real(0.0);

    const auto nlevels = static_cast<int>(geom.size());

#ifdef AMREX_FEATURE_HYPRE_SSAMG
    auto hypre_solver_id = use_hypre_ssamg ? HypreSolverID::SSAMG
                                           : HypreSolverID::BoomerAMG;
#else
    auto hypre_solver_id = HypreSolverID::BoomerAMG;
#endif

    if (prob_type == 1) { // Poisson
        if (composite_solve) {
            HypreMLABecLap hypre_mlabeclap(geom, grids, dmap, hypre_solver_id);
            hypre_mlabeclap.setVerbose(verbose);

            hypre_mlabeclap.setup(Real(0.0), Real(-1.0), {}, {},
                                  {AMREX_D_DECL(LinOpBCType::Dirichlet,
                                                LinOpBCType::Dirichlet,
                                                LinOpBCType::Dirichlet)},
                                  {AMREX_D_DECL(LinOpBCType::Dirichlet,
                                                LinOpBCType::Dirichlet,
                                                LinOpBCType::Dirichlet)},
                                  GetVecOfConstPtrs(solution));

            hypre_mlabeclap.solve(GetVecOfPtrs(solution), GetVecOfConstPtrs(rhs),
                                  tol_rel, tol_abs);
        } else {
            for (int ilev = 0; ilev < nlevels; ++ilev) {
                HypreMLABecLap hypre_mlabeclap({geom[ilev]}, {grids[ilev]}, {dmap[ilev]}, hypre_solver_id);
                hypre_mlabeclap.setVerbose(verbose);

                std::pair<MultiFab const*, IntVect> coarse_bc{nullptr,IntVect(0)};
                if (ilev > 0) {
                    coarse_bc.first = &solution[ilev-1];
                    coarse_bc.second = IntVect(ref_ratio);
                }

                hypre_mlabeclap.setup(Real(0.0), Real(-1.0), {}, {},
                                      {AMREX_D_DECL(LinOpBCType::Dirichlet,
                                                    LinOpBCType::Dirichlet,
                                                    LinOpBCType::Dirichlet)},
                                      {AMREX_D_DECL(LinOpBCType::Dirichlet,
                                                    LinOpBCType::Dirichlet,
                                                    LinOpBCType::Dirichlet)},
                                      {&solution[ilev]},
                                      coarse_bc);

                hypre_mlabeclap.solve({&solution[ilev]}, {&rhs[ilev]}, tol_rel, tol_abs);
            }
        }
    } else if (prob_type == 2) { // ABecLaplacian
        Vector<Array<MultiFab,AMREX_SPACEDIM>> face_bcoef(nlevels);
        for (int ilev = 0; ilev < nlevels; ++ilev) {
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                const BoxArray& ba = amrex::convert(bcoef[ilev].boxArray(),
                                                    IntVect::TheDimensionVector(idim));
                face_bcoef[ilev][idim].define(ba, bcoef[ilev].DistributionMap(), 1, 0);
            }
            amrex::average_cellcenter_to_face(GetArrOfPtrs(face_bcoef[ilev]),
                                              bcoef[ilev], geom[ilev]);
        }

        if (composite_solve) {
            HypreMLABecLap hypre_mlabeclap(geom, grids, dmap, hypre_solver_id);
            hypre_mlabeclap.setVerbose(verbose);

            hypre_mlabeclap.setup(ascalar, bscalar,
                                  GetVecOfConstPtrs(acoef),
                                  GetVecOfArrOfConstPtrs(face_bcoef),
                                  {AMREX_D_DECL(LinOpBCType::Dirichlet,
                                                LinOpBCType::Neumann,
                                                LinOpBCType::Neumann)},
                                  {AMREX_D_DECL(LinOpBCType::Neumann,
                                                LinOpBCType::Dirichlet,
                                                LinOpBCType::Neumann)},
                                  GetVecOfConstPtrs(solution));

            hypre_mlabeclap.solve(GetVecOfPtrs(solution), GetVecOfConstPtrs(rhs),
                                  tol_rel, tol_abs);
        } else {
            for (int ilev = 0; ilev < nlevels; ++ilev) {
                HypreMLABecLap hypre_mlabeclap({geom[ilev]}, {grids[ilev]}, {dmap[ilev]}, hypre_solver_id);
                hypre_mlabeclap.setVerbose(verbose);

                std::pair<MultiFab const*, IntVect> coarse_bc{nullptr,IntVect(0)};
                if (ilev > 0) {
                    coarse_bc.first = &solution[ilev-1];
                    coarse_bc.second = IntVect(ref_ratio);
                }

                hypre_mlabeclap.setup(ascalar, bscalar,
                                      {&acoef[ilev]},
                                      {GetArrOfConstPtrs(face_bcoef[ilev])},
                                      {AMREX_D_DECL(LinOpBCType::Dirichlet,
                                                    LinOpBCType::Neumann,
                                                    LinOpBCType::Neumann)},
                                      {AMREX_D_DECL(LinOpBCType::Neumann,
                                                    LinOpBCType::Dirichlet,
                                                    LinOpBCType::Neumann)},
                                      {&solution[ilev]},
                                      coarse_bc);

                hypre_mlabeclap.solve({&solution[ilev]}, {&rhs[ilev]}, tol_rel, tol_abs);
            }
        }
    } else {
        amrex::Abort("Unsupported prob_type: " + std::to_string(prob_type));
    }
}
#endif
