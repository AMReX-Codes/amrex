
#include <AMReX_HypreABecLap.H>
#include <string>
#include <algorithm>

#include <AMReX_Habec_K.H>

#include <_hypre_struct_mv.h>

namespace amrex {

HypreABecLap::HypreABecLap(const BoxArray& grids, const DistributionMapping& dmap,
                           const Geometry& geom_, MPI_Comm comm_)
    : Hypre(grids, dmap, geom_, comm_)
{
}

HypreABecLap::~HypreABecLap ()
{
    HYPRE_StructPFMGDestroy(solver);
    solver = NULL;
    HYPRE_StructMatrixDestroy(A);
    A = NULL;
//    HYPRE_StructVectorDestroy(b);
//    b = NULL;
//    HYPRE_StructVectorDestroy(x);
//    x = NULL;
    HYPRE_StructGridDestroy(grid);
    grid = NULL;
}

void
HypreABecLap::solve(MultiFab& soln, const MultiFab& rhs, Real reltol, Real abstol,
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

    // do this repeatedly to avoid memory leak
    HYPRE_StructVectorCreate(comm, grid, &b);
    HYPRE_StructVectorCreate(comm, grid, &x);
    //
    HYPRE_StructVectorInitialize(b);
    HYPRE_StructVectorInitialize(x);
    //
    loadVectors(soln, rhs);
    //
    HYPRE_StructVectorAssemble(x); 
    HYPRE_StructVectorAssemble(b); 

    HYPRE_StructPFMGSetMaxIter(solver, maxiter);
    HYPRE_StructPFMGSetTol(solver, reltol);
    if (abstol > 0.0)
    {
        Real bnorm;
        bnorm = hypre_StructInnerProd((hypre_StructVector *) b,
                                      (hypre_StructVector *) b);
        bnorm = std::sqrt(bnorm);

        const BoxArray& grids = rhs.boxArray();
        Real volume = grids.d_numPts();
        Real reltol_new = (bnorm > 0.0) ? (abstol / bnorm * std::sqrt(volume)) : reltol;

        if (reltol_new > reltol) {
            HYPRE_StructPFMGSetTol(solver, reltol_new);
        }
    }

    HYPRE_StructPFMGSolve(solver, A, b, x);

    if (verbose >= 2)
    {
        HYPRE_Int num_iterations;
        Real res;
        HYPRE_StructPFMGGetNumIterations(solver, &num_iterations);
        HYPRE_StructPFMGGetFinalRelativeResidualNorm(solver, &res);
        
        amrex::Print() << "\n" << num_iterations
                       << " Hypre PFMG Iterations, Relative Residual "
                       << res << std::endl;
    }

    getSolution(soln);

    // do this repeatedly to avoid memory leak
    HYPRE_StructVectorDestroy(b);
    b = NULL;
    HYPRE_StructVectorDestroy(x);
    x = NULL;
}

void
HypreABecLap::getSolution (MultiFab& soln)
{
    FArrayBox rfab;
    for (MFIter mfi(soln); mfi.isValid(); ++mfi)
    {
        const Box &reg = mfi.validbox();

        FArrayBox *xfab;
        if (soln.nGrow() == 0) { // need a temporary if soln is the wrong size
            xfab = &soln[mfi];
        }
        else {
            xfab = &rfab;
            xfab->resize(reg);
        }
        
        auto reglo = Hypre::loV(reg);
        auto reghi = Hypre::hiV(reg);
        HYPRE_StructVectorGetBoxValues(x, reglo.data(), reghi.data(), xfab->dataPtr());

        if (soln.nGrow() != 0) {
            soln[mfi].copy<RunOn::Host>(*xfab, 0, 0, 1);
        }
    }
}

void
HypreABecLap::prepareSolver ()
{
    BL_PROFILE("HypreABecLap::prepareSolver()");

    const BoxArray& ba = acoefs.boxArray();
    const DistributionMapping& dm = acoefs.DistributionMap();

    HYPRE_StructGridCreate(comm, AMREX_SPACEDIM, &grid);

    Array<HYPRE_Int,AMREX_SPACEDIM> is_periodic {AMREX_D_DECL(0,0,0)};
    for (int i = 0; i < AMREX_SPACEDIM; i++) {
        if (geom.isPeriodic(i)) {
            is_periodic[i] = geom.period(i);
            AMREX_ASSERT(Hypre::ispow2(is_periodic[i]));
            AMREX_ASSERT(geom.Domain().smallEnd(i) == 0);
        }
    }
    if (geom.isAnyPeriodic()) {
        HYPRE_StructGridSetPeriodic(grid, is_periodic.data());
    }

    for (MFIter mfi(acoefs); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();
        auto lo = Hypre::loV(bx);
        auto hi = Hypre::hiV(bx);
        HYPRE_StructGridSetExtents(grid, lo.data(), hi.data());
    }

    HYPRE_StructGridAssemble(grid);

#if (AMREX_SPACEDIM == 2)
    HYPRE_Int offsets[regular_stencil_size][2] = {{ 0,  0},    // 0
                                                  {-1,  0},    // 1
                                                  { 1,  0},    // 2
                                                  { 0, -1},    // 3
                                                  { 0,  1}};   // 4
#elif (AMREX_SPACEDIM == 3)
    HYPRE_Int offsets[regular_stencil_size][3] = {{ 0,  0,  0},   // 0
                                                  {-1,  0,  0},   // 1
                                                  { 1,  0,  0},   // 2
                                                  { 0, -1,  0},   // 3
                                                  { 0,  1,  0},   // 4
                                                  { 0,  0, -1},   // 5
                                                  { 0,  0,  1}};  // 6
#endif

    HYPRE_StructStencil stencil;
    HYPRE_StructStencilCreate(AMREX_SPACEDIM, regular_stencil_size, &stencil);

    for (int i = 0; i < regular_stencil_size; ++i) {
        HYPRE_StructStencilSetElement(stencil, i, offsets[i]);
    }

    HYPRE_StructMatrixCreate(comm, grid, stencil, &A);
    HYPRE_StructMatrixInitialize(A);

    HYPRE_StructStencilDestroy(stencil); 

    HYPRE_StructVectorCreate(comm, grid, &b);
    HYPRE_StructVectorCreate(comm, grid, &x);

    HYPRE_StructVectorInitialize(b);
    HYPRE_StructVectorInitialize(x);

    // A.SetValues() & A.assemble()
    Array<HYPRE_Int,regular_stencil_size> stencil_indices;
    std::iota(stencil_indices.begin(), stencil_indices.end(), 0);
    const Real* dx = geom.CellSize();
    const int bho = (m_maxorder > 2) ? 1 : 0;
    BaseFab<GpuArray<Real,regular_stencil_size>> rfab;
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

        // add b.c.'s for A matrix
        for (OrientationIter oit; oit; oit++)
        {
            Orientation ori = oit();
            int cdir(ori);
            int idim = ori.coordDir();
            const int bctype = bcs_i[cdir][0];
            const Real  &bcl = bcl_i[cdir];
            const Mask  &msk = m_bndry->bndryMasks(ori)[mfi];
   
            amrex_hpmat(reg, rfab, bcoefs[idim][mfi], msk, scalar_b, dx, cdir, bctype, bcl, bho);   
        }
        
        amrex_hpdiag(reg, rfab, diaginv[mfi]);
        Real* mat = (Real*) rfab.dataPtr();

        auto reglo = Hypre::loV(reg);
        auto reghi = Hypre::hiV(reg);
        HYPRE_StructMatrixSetBoxValues(A, reglo.data(), reghi.data(),
                                       regular_stencil_size, stencil_indices.data(),
                                       mat);
    }
    HYPRE_StructMatrixAssemble(A);

    // create solver
    HYPRE_StructPFMGCreate(comm, &solver);
    int logging = (verbose >= 2) ? 1 : 0;
    HYPRE_StructPFMGSetLogging(solver, logging);
    HYPRE_StructPFMGSetup(solver, A, b, x);

    HYPRE_StructVectorDestroy(b);
    b = NULL;
    HYPRE_StructVectorDestroy(x);
    x = NULL;
}


void
HypreABecLap::loadVectors (MultiFab& soln, const MultiFab& rhs)
{
    BL_PROFILE("HypreABecLap::loadVectors()");

    soln.setVal(0.0);

    FArrayBox rhsfab;
    for (MFIter mfi(soln); mfi.isValid(); ++mfi)
    {
        const Box &reg = mfi.validbox();

        auto reglo = Hypre::loV(reg);
        auto reghi = Hypre::hiV(reg);
        HYPRE_StructVectorSetBoxValues(x, reglo.data(), reghi.data(), soln[mfi].dataPtr());
        
        rhsfab.resize(reg);
        rhsfab.copy<RunOn::Host>(rhs[mfi],reg);
        rhsfab.mult<RunOn::Host>(diaginv[mfi]);

        HYPRE_StructVectorSetBoxValues(b, reglo.data(), reghi.data(), rhsfab.dataPtr());
    }
}

}
