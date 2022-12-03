#ifndef MY_TEST_H_
#define MY_TEST_H_

#include <AMReX_MLMG.H>
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLPoisson.H>
#include <AMReX_MultiFabUtil.H>

class MyTest
{
public:

    MyTest ();

    void solve ();
    void writePlotfile () const;

public: // make these public for cuda
    void initProbPoisson ();
    void initProbABecLaplacian ();

private:

    void readParameters ();
    void initData ();

    template <typename MF>
    void solvePoisson ();

    template <typename MF>
    void solveABecLaplacian ();

    int max_level = 1;
    int ref_ratio = 2;
    int n_cell = 128;
    int max_grid_size = 64;

    bool composite_solve = false;

    int prob_type = 1;  // 1. Poisson,  2. ABecLaplacian

    bool single_precision = true;

    // For MLMG solver
    int verbose = 2;
    int bottom_verbose = 0;
    int max_iter = 100;
    int max_fmg_iter = 0;
    int linop_maxorder = 2;
    bool agglomeration = true;
    bool consolidation = true;
    int max_coarsening_level = 30;

    amrex::Vector<amrex::Geometry> geom;
    amrex::Vector<amrex::BoxArray> grids;
    amrex::Vector<amrex::DistributionMapping> dmap;

    amrex::Vector<amrex::MultiFab> solution;
    amrex::Vector<amrex::MultiFab> rhs;
    amrex::Vector<amrex::MultiFab> exact_solution;
    amrex::Vector<amrex::MultiFab> acoef;
    amrex::Vector<amrex::MultiFab> bcoef;

    amrex::Real ascalar = 1.e-3;
    amrex::Real bscalar = 1.0;
};

template <typename MF>
void
MyTest::solvePoisson ()
{
    using namespace amrex;

    LPInfo info;
    info.setAgglomeration(agglomeration);
    info.setConsolidation(consolidation);
    info.setMaxCoarseningLevel(max_coarsening_level);

    using T = typename MF::value_type;
    const T tol_rel = std::is_same<double,typename MF::value_type>::value ?
        T(1.e-10) : T(1.e-4);
    const Real tol_abs = 0.0;

    const int nlevels = geom.size();

    if (composite_solve)
    {

        MLPoissonT<MF> mlpoisson(geom, grids, dmap, info);

        mlpoisson.setMaxOrder(linop_maxorder);

        // This is a problem with Dirichlet BC
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

        MLMGT<MF> mlmg(mlpoisson);
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
            MLPoissonT<MF> mlpoisson({geom[ilev]}, {grids[ilev]}, {dmap[ilev]}, info);

            mlpoisson.setMaxOrder(linop_maxorder);

            // This is a problem with Dirichlet BC
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

            MLMGT<MF> mlmg(mlpoisson);
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

template <typename MF>
void
MyTest::solveABecLaplacian ()
{
    using namespace amrex;

    LPInfo info;
    info.setAgglomeration(agglomeration);
    info.setConsolidation(consolidation);
    info.setMaxCoarseningLevel(max_coarsening_level);

    using T = typename MF::value_type;
    const T tol_rel = std::is_same<double,typename MF::value_type>::value ?
        T(1.e-10) : T(1.e-4);
    const T tol_abs = T(0.0);

    const int nlevels = geom.size();

    if (composite_solve)
    {
        MLABecLaplacianT<MF> mlabec(geom, grids, dmap, info);

        mlabec.setMaxOrder(linop_maxorder);

        // This is a problem with Dirichlet BC
        mlabec.setDomainBC({AMREX_D_DECL(LinOpBCType::Dirichlet,
                                         LinOpBCType::Dirichlet,
                                         LinOpBCType::Dirichlet)},
                           {AMREX_D_DECL(LinOpBCType::Dirichlet,
                                         LinOpBCType::Dirichlet,
                                         LinOpBCType::Dirichlet)});

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

        MLMGT<MF> mlmg(mlabec);
        mlmg.setMaxIter(max_iter);
        mlmg.setMaxFmgIter(max_fmg_iter);
        mlmg.setVerbose(verbose);
        mlmg.setBottomVerbose(bottom_verbose);
        mlmg.solve(GetVecOfPtrs(solution), GetVecOfConstPtrs(rhs), tol_rel, tol_abs);
    }
    else
    {
        for (int ilev = 0; ilev < nlevels; ++ilev)
        {
            MLABecLaplacianT<MF> mlabec({geom[ilev]}, {grids[ilev]}, {dmap[ilev]}, info);

            mlabec.setMaxOrder(linop_maxorder);

            // This is a 3d problem with Dirichlet BC
            mlabec.setDomainBC({AMREX_D_DECL(LinOpBCType::Dirichlet,
                                             LinOpBCType::Dirichlet,
                                             LinOpBCType::Dirichlet)},
                               {AMREX_D_DECL(LinOpBCType::Dirichlet,
                                             LinOpBCType::Dirichlet,
                                             LinOpBCType::Dirichlet)});

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

            MLMGT<MF> mlmg(mlabec);
            mlmg.setMaxIter(max_iter);
            mlmg.setMaxFmgIter(max_fmg_iter);
            mlmg.setVerbose(verbose);
            mlmg.setBottomVerbose(bottom_verbose);
            mlmg.solve({&solution[ilev]}, {&rhs[ilev]}, tol_rel, tol_abs);
        }
    }
}

#endif
