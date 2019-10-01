#include "MyTest.H"

#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLPoisson.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>

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
        solvePoisson();
    } else if (prob_type == 2) {
        solveABecLaplacian();
    } else if (prob_type == 3) {
        solveABecLaplacianInhomNeumann();
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

    const Real tol_rel = 1.e-10;
    const Real tol_abs = 0.0;

    const int nlevels = geom.size();

    if (composite_solve)
    {

        MLPoisson mlpoisson(geom, grids, dmap, info);

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
    info.setMaxCoarseningLevel(max_coarsening_level);

    const Real tol_rel = 1.e-10;
    const Real tol_abs = 0.0;

    const int nlevels = geom.size();

    if (composite_solve)
    {

        MLABecLaplacian mlabec(geom, grids, dmap, info);

        mlabec.setMaxOrder(linop_maxorder);

        // This is a 3d problem with homogeneous Neumann BC
        mlabec.setDomainBC({AMREX_D_DECL(LinOpBCType::Neumann,
                                         LinOpBCType::Neumann,
                                         LinOpBCType::Neumann)},
                           {AMREX_D_DECL(LinOpBCType::Neumann,
                                         LinOpBCType::Neumann,
                                         LinOpBCType::Neumann)});

        for (int ilev = 0; ilev < nlevels; ++ilev)
        {
            // for problem with pure homogeneous Neumann BC, we could pass a nullptr
            mlabec.setLevelBC(ilev, nullptr);
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
            
            mlabec.setMaxOrder(linop_maxorder);
            
            // This is a 3d problem with homogeneous Neumann BC
            mlabec.setDomainBC({AMREX_D_DECL(LinOpBCType::Neumann,
                                             LinOpBCType::Neumann,
                                             LinOpBCType::Neumann)},
                               {AMREX_D_DECL(LinOpBCType::Neumann,
                                             LinOpBCType::Neumann,
                                             LinOpBCType::Neumann)});
            
            if (ilev > 0) {
                mlabec.setCoarseFineBC(&solution[ilev-1], ref_ratio);
            }

            // for problem with pure homogeneous Neumann BC, we could pass a nullptr
            mlabec.setLevelBC(0, nullptr);

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

    // Since this problem has Neumann BC, solution + constant is also a
    // solution.  So we are going to shift the solution by a constant
    // for comparison with the "exact solution".
    const Real npts = grids[0].d_numPts();
    const Real avg1 = exact_solution[0].sum();
    const Real avg2 = solution[0].sum();
    const Real offset = (avg1-avg2)/npts;
    for (int ilev = 0; ilev < nlevels; ++ilev) {
        solution[ilev].plus(offset, 0, 1, 0);
    }
}

void
MyTest::solveABecLaplacianInhomNeumann ()
{
    LPInfo info;
    info.setAgglomeration(agglomeration);
    info.setConsolidation(consolidation);
    info.setMaxCoarseningLevel(max_coarsening_level);

    const Real tol_rel = 1.e-10;
    const Real tol_abs = 0.0;

    const int nlevels = geom.size();

    if (composite_solve)
    {

        MLABecLaplacian mlabec(geom, grids, dmap, info);

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

    // Since this problem has Neumann BC, solution + constant is also a
    // solution.  So we are going to shift the solution by a constant
    // for comparison with the "exact solution".
    const Real npts = grids[0].d_numPts();
    const Real avg1 = exact_solution[0].sum();
    const Real avg2 = solution[0].sum();
    const Real offset = (avg1-avg2)/npts;
    for (int ilev = 0; ilev < nlevels; ++ilev) {
        solution[ilev].plus(offset, 0, 1, 0);
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

    pp.query("verbose", verbose);
    pp.query("bottom_verbose", bottom_verbose);
    pp.query("max_iter", max_iter);
    pp.query("max_fmg_iter", max_fmg_iter);
    pp.query("linop_maxorder", linop_maxorder);
    pp.query("agglomeration", agglomeration);
    pp.query("consolidation", consolidation);
    pp.query("max_coarsening_level", max_coarsening_level);

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
    
    if (prob_type == 2 or prob_type == 3) {
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
    } else if (prob_type == 3) {
        initProbABecLaplacianInhomNeumann();
    } else {
        amrex::Abort("Unknown prob_type "+std::to_string(prob_type));
    }
}

