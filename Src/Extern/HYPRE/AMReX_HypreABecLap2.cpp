#include <AMReX_HypreABecLap2.H>

#include <AMReX_Habec_K.H>

#include <cmath>
#include <numeric>
#include <algorithm>
#include <type_traits>

namespace amrex
{

HypreABecLap2::HypreABecLap2 (const BoxArray& grids, const DistributionMapping& dmap,
                              const Geometry& geom_, MPI_Comm comm_)
    : Hypre(grids, dmap, geom_, comm_)
{
}

HypreABecLap2::~HypreABecLap2 ()
{
    HYPRE_BoomerAMGDestroy(solver);
    solver = NULL;
    HYPRE_SStructMatrixDestroy(A);
    A = NULL;
//    HYPRE_SStructVectorDestroy(b);  // done in solve()
//    b = NULL;
//    HYPRE_SStructVectorDestroy(x);
//    x = NULL;
    HYPRE_SStructGraphDestroy(graph);
    graph = NULL;
    HYPRE_SStructStencilDestroy(stencil);
    stencil = NULL;
    HYPRE_SStructGridDestroy(hgrid);
    hgrid = NULL;
}

void
HypreABecLap2::solve (MultiFab& soln, const MultiFab& rhs, Real reltol, Real abstol, 
                      int maxiter, const BndryData& bndry, int max_bndry_order)
{
    if (solver == NULL || m_bndry != &bndry || m_maxorder != max_bndry_order)
    {
        m_bndry = &bndry;
        m_maxorder = max_bndry_order;
        m_factory = &(rhs.Factory());
        prepareSolver();
    }
    else
    {
        m_factory = &(rhs.Factory());
    }

    // We have to do this repeatedly to avoid memory leak due to Hypre bug
    HYPRE_SStructVectorCreate(comm, hgrid, &b);
    HYPRE_SStructVectorSetObjectType(b, HYPRE_PARCSR);
    HYPRE_SStructVectorInitialize(b);
    //
    HYPRE_SStructVectorCreate(comm, hgrid, &x);
    HYPRE_SStructVectorSetObjectType(x, HYPRE_PARCSR);
    HYPRE_SStructVectorInitialize(x);
    //
    loadVectors(soln, rhs);
    //
    HYPRE_SStructVectorAssemble(b);
    HYPRE_SStructVectorAssemble(x);

    HYPRE_BoomerAMGSetMinIter(solver, 1);
    HYPRE_BoomerAMGSetMaxIter(solver, maxiter);
    HYPRE_BoomerAMGSetTol(solver, reltol);
    if (abstol > 0.0)
    {
        Real bnorm;
        hypre_SStructInnerProd((hypre_SStructVector *) b,
                               (hypre_SStructVector *) b,
                               &bnorm);
        bnorm = std::sqrt(bnorm);

        const BoxArray& grids = acoefs.boxArray();
        Real volume = grids.d_numPts();
        Real reltol_new = bnorm > 0.0 ? abstol / (bnorm+1.e-100) * std::sqrt(volume) : reltol;

        if (reltol_new > reltol) {
            HYPRE_BoomerAMGSetTol(solver, reltol_new);
        }
    }

    HYPRE_ParCSRMatrix par_A;
    HYPRE_ParVector par_b;
    HYPRE_ParVector par_x;
    
    HYPRE_SStructMatrixGetObject(A, (void**) &par_A);
    HYPRE_SStructVectorGetObject(b, (void**) &par_b);
    HYPRE_SStructVectorGetObject(x, (void**) &par_x);

    HYPRE_BoomerAMGSolve(solver, par_A, par_b, par_x);

    if (verbose >= 2)
    {
        HYPRE_Int num_iterations;
        Real res;
        HYPRE_BoomerAMGGetNumIterations(solver, &num_iterations);
        HYPRE_BoomerAMGGetFinalRelativeResidualNorm(solver, &res);

        amrex::Print() << "\n" << num_iterations
                       << " Hypre SS BoomerAMG Iterations, Relative Residual "
                       << res << std::endl;
    }

    getSolution(soln);

    // We have to do this repeatedly to avoid memory leak due to Hypre bug
    HYPRE_SStructVectorDestroy(b);
    b = NULL;
    HYPRE_SStructVectorDestroy(x);
    x = NULL;
}

void
HypreABecLap2::getSolution (MultiFab& soln)
{
    HYPRE_SStructVectorGather(x);

    const HYPRE_Int part = 0;

    FArrayBox fab;
    for (MFIter mfi(soln); mfi.isValid(); ++mfi)
    {
        const Box &reg = mfi.validbox();

        FArrayBox *xfab;
        if (soln.nGrow() == 0) { // need a temporary if soln is the wrong size
            xfab = &soln[mfi];
        }
        else {
            xfab = &fab;
            xfab->resize(reg);
        }

        auto reglo = Hypre::loV(reg);
        auto reghi = Hypre::hiV(reg);
        HYPRE_SStructVectorGetBoxValues(x, part, reglo.data(), reghi.data(),
                                        0, xfab->dataPtr());

        if (soln.nGrow() != 0) {
            soln[mfi].copy<RunOn::Host>(*xfab, 0, 0, 1);
        }
    }
}

void
HypreABecLap2::prepareSolver ()
{
    BL_PROFILE("HypreABecLap2::prepareSolver()");

    HYPRE_SStructGridCreate(comm, AMREX_SPACEDIM, 1, &hgrid);

    Array<HYPRE_Int,AMREX_SPACEDIM> is_periodic {AMREX_D_DECL(0,0,0)};
    for (int i = 0; i < AMREX_SPACEDIM; i++) {
        if (geom.isPeriodic(i)) {
            is_periodic[i] = geom.period(i);
            AMREX_ASSERT(Hypre::ispow2(is_periodic[i]));
            AMREX_ASSERT(geom.Domain().smallEnd(i) == 0);
        }
    }
    if (geom.isAnyPeriodic()) {
        HYPRE_SStructGridSetPeriodic(hgrid, 0, is_periodic.data());
    }

    for (MFIter mfi(acoefs); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();
        auto lo = Hypre::loV(bx);
        auto hi = Hypre::hiV(bx);
        HYPRE_SStructGridSetExtents(hgrid, 0, lo.data(), hi.data());
    }

    // All variables are cell-centered
    HYPRE_SStructVariable vars[1] = {HYPRE_SSTRUCT_VARIABLE_CELL};
    HYPRE_SStructGridSetVariables(hgrid, 0, 1, vars);

    HYPRE_SStructGridAssemble(hgrid);

    // Setup stencils
#if (AMREX_SPACEDIM == 2)
    HYPRE_Int offsets[regular_stencil_size][2] = {{ 0,  0},
                                                  {-1,  0},
                                                  { 1,  0},
                                                  { 0, -1},
                                                  { 0,  1}};
#elif (AMREX_SPACEDIM == 3)
    HYPRE_Int offsets[regular_stencil_size][3] = {{ 0,  0,  0},
                                                  {-1,  0,  0},
                                                  { 1,  0,  0},
                                                  { 0, -1,  0},
                                                  { 0,  1,  0},
                                                  { 0,  0, -1},
                                                  { 0,  0,  1}};
#endif

    HYPRE_SStructStencilCreate(AMREX_SPACEDIM, regular_stencil_size, &stencil);

    for (int i = 0; i < regular_stencil_size; i++) {
        HYPRE_SStructStencilSetEntry(stencil, i, offsets[i], 0);
    }

    HYPRE_SStructGraphCreate(comm, hgrid, &graph);
    HYPRE_SStructGraphSetObjectType(graph, HYPRE_PARCSR);

    HYPRE_SStructGraphSetStencil(graph, 0, 0, stencil);

    HYPRE_SStructGraphAssemble(graph);

    HYPRE_SStructMatrixCreate(comm, graph, &A);
    HYPRE_SStructMatrixSetObjectType(A, HYPRE_PARCSR);
    HYPRE_SStructMatrixInitialize(A);

    // A.SetValues() & A.assemble()
    Array<HYPRE_Int,regular_stencil_size> stencil_indices;
    std::iota(stencil_indices.begin(), stencil_indices.end(), 0);
    const HYPRE_Int part = 0;
    const Real* dx = geom.CellSize();
    const int bho = (m_maxorder > 2) ? 1 : 0;
    BaseFab<GpuArray<Real, regular_stencil_size>> rfab;
    for (MFIter mfi(acoefs); mfi.isValid(); ++mfi)
    {
        const Box &reg = mfi.validbox();

        rfab.resize(reg);
        amrex_hpacoef(reg, rfab, acoefs[mfi], scalar_a);
         
        for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
            amrex_hpbcoef(reg, rfab, bcoefs[idim][mfi], scalar_b, dx, idim);
        }

        const Vector< Vector<BoundCond> > & bcs_i = m_bndry->bndryConds(mfi);
        const BndryData::RealTuple        & bcl_i = m_bndry->bndryLocs(mfi);
        
        for (OrientationIter oit; oit; oit++)
        {
            Orientation ori = oit();
            int cdir(ori);
            int idim = ori.coordDir();
            const int bctype = bcs_i[cdir][0];
            const Real &bcl  = bcl_i[cdir];
            const Mask &msk  = m_bndry->bndryMasks(ori)[mfi];

            amrex_hpmat(reg, rfab, bcoefs[idim][mfi], msk, scalar_b, dx, cdir, bctype, bcl, bho);
        }

        amrex_hpdiag(reg, rfab, diaginv[mfi]); 
        Real* mat = (Real*) rfab.dataPtr();

        // initialize matrix
        auto reglo = Hypre::loV(reg);
        auto reghi = Hypre::hiV(reg);
        HYPRE_SStructMatrixSetBoxValues(A, part, reglo.data(), reghi.data(),
                                        0, regular_stencil_size, stencil_indices.data(),
                                        mat);
    }    
    HYPRE_SStructMatrixAssemble(A);   

    // create solver
    HYPRE_BoomerAMGCreate(&solver);

    HYPRE_BoomerAMGSetOldDefault(solver); // Falgout coarsening with modified classical interpolation
//    HYPRE_BoomerAMGSetCoarsenType(solver, 6);
//    HYPRE_BoomerAMGSetCycleType(solver, 1);
    HYPRE_BoomerAMGSetRelaxType(solver, 6);   /* G-S/Jacobi hybrid relaxation */
    HYPRE_BoomerAMGSetRelaxOrder(solver, 1);   /* uses C/F relaxation */
    HYPRE_BoomerAMGSetNumSweeps(solver, 2);   /* Sweeeps on each level */
//    HYPRE_BoomerAMGSetStrongThreshold(solver, 0.6); // default is 0.25

    int logging = (verbose >= 2) ? 1 : 0;
    HYPRE_BoomerAMGSetLogging(solver, logging);

    HYPRE_ParCSRMatrix par_A;
    HYPRE_SStructMatrixGetObject(A, (void**) &par_A);
    HYPRE_BoomerAMGSetup(solver, par_A, NULL, NULL);
}

void
HypreABecLap2::loadVectors (MultiFab& soln, const MultiFab& rhs)
{
    BL_PROFILE("HypreABecLap2::loadVectors()");

    soln.setVal(0.0);

    const HYPRE_Int part = 0;
    FArrayBox rhsfab;
    for (MFIter mfi(soln); mfi.isValid(); ++mfi)
    {
        const Box &reg = mfi.validbox();

        auto reglo = Hypre::loV(reg);
        auto reghi = Hypre::hiV(reg);
        HYPRE_SStructVectorSetBoxValues(x, part, reglo.data(), reghi.data(),
                                        0, soln[mfi].dataPtr());

        rhsfab.resize(reg);
        rhsfab.copy<RunOn::Host>(rhs[mfi],reg);
        rhsfab.mult<RunOn::Host>(diaginv[mfi]);

        HYPRE_SStructVectorSetBoxValues(b, part, reglo.data(), reghi.data(),
                                        0, rhsfab.dataPtr());
    }
}

}
