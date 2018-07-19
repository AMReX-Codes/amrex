
#include <AMReX_Hypre.H>
#include <AMReX_HypreABec_F.H>
#include <string>
#include <algorithm>

#include <_hypre_struct_mv.h>

namespace amrex {

constexpr int HypreABecLap::stencil_size;
    
HypreABecLap::HypreABecLap(const BoxArray& grids,
                           const DistributionMapping& dmap,
                           const Geometry& geom_,
                           MPI_Comm comm_)
    : comm(comm_),
      verbose(0),
      geom(geom_)
{
    int num_procs, myid;
    MPI_Comm_size(comm, &num_procs );
    MPI_Comm_rank(comm, &myid );

    HYPRE_StructGridCreate(comm, AMREX_SPACEDIM, &grid);

    for (int i = 0; i < AMREX_SPACEDIM; i++) {
        is_periodic[i] = 0;
        if (geom.isPeriodic(i)) {
            is_periodic[i] = geom.period(i);
            BL_ASSERT(Hypre::ispow2(is_periodic[i]));
            BL_ASSERT(geom.Domain().smallEnd(i) == 0);
        }
    }
    if (geom.isAnyPeriodic()) {
        HYPRE_StructGridSetPeriodic(grid, is_periodic);
    }

    if (num_procs != 1) {
        for (int i = 0; i < grids.size(); i++) {
            if (dmap[i] == myid) {
                HYPRE_StructGridSetExtents(grid, Hypre::loV(grids[i]), Hypre::hiV(grids[i]));
            }
        }
    }
    else {
        for (int i = 0; i < grids.size(); i++) {
            HYPRE_StructGridSetExtents(grid, Hypre::loV(grids[i]), Hypre::hiV(grids[i]));
        }
    }

    HYPRE_StructGridAssemble(grid);

#if (AMREX_SPACEDIM == 2)
    int offsets[stencil_size][2] = {{ 0,  0},    // 0
                                    {-1,  0},    // 1
                                    { 1,  0},    // 2
                                    { 0, -1},    // 3
                                    { 0,  1}};   // 4
#elif (AMREX_SPACEDIM == 3)
    int offsets[stencil_size][3] = {{ 0,  0,  0},   // 0
                                    {-1,  0,  0},   // 1
                                    { 1,  0,  0},   // 2
                                    { 0, -1,  0},   // 3
                                    { 0,  1,  0},   // 4
                                    { 0,  0, -1},   // 5
                                    { 0,  0,  1}};  // 6
#endif

    HYPRE_StructStencil stencil;

    HYPRE_StructStencilCreate(AMREX_SPACEDIM, stencil_size, &stencil);

    for (int i = 0; i < stencil_size; i++) {
        HYPRE_StructStencilSetElement(stencil, i, offsets[i]);
    }

    HYPRE_StructMatrixCreate(comm, grid, stencil, &A);
    HYPRE_StructMatrixInitialize(A);

    HYPRE_StructStencilDestroy(stencil); 

    HYPRE_StructVectorCreate(comm, grid, &b);
    HYPRE_StructVectorCreate(comm, grid, &x);

    HYPRE_StructVectorInitialize(b);
    HYPRE_StructVectorInitialize(x);

    int ncomp=1;
    int ngrow=0;
    acoefs.reset(new MultiFab(grids, dmap, ncomp, ngrow));
    acoefs->setVal(0.0);
 
    for (int i = 0; i < AMREX_SPACEDIM; i++) {
        BoxArray edge_boxes(grids);
        edge_boxes.surroundingNodes(i);
        bcoefs[i].reset(new MultiFab(edge_boxes, dmap, ncomp, ngrow));
    }
}

HypreABecLap::~HypreABecLap ()
{
    HYPRE_StructGridDestroy(grid);
    HYPRE_StructMatrixDestroy(A);
    HYPRE_StructVectorDestroy(b);
    HYPRE_StructVectorDestroy(x);
}

void
HypreABecLap::setScalars(Real sa, Real sb)
{
    scalar_a = sa;
    scalar_b = sb;
}

void
HypreABecLap::setACoeffs(const MultiFab& alpha)
{
    MultiFab::Copy(*acoefs, alpha, 0, 0, 1, 0);
}

void
HypreABecLap::setBCoeffs(const std::array<const MultiFab*,AMREX_SPACEDIM>& beta)
{
    for (int idim=0; idim<AMREX_SPACEDIM; idim++) {
        MultiFab::Copy(*bcoefs[idim], *beta[idim], 0, 0, 1, 0);
    }
}

void
HypreABecLap::setVerbose(int _verbose)
{
  verbose = _verbose;
}

void
HypreABecLap::solve(MultiFab& soln, const MultiFab& rhs, Real reltol, Real abstol,
                    int maxiter, const BndryData& bndry, int max_bndry_order)
{
    // set up solver
    const BoxArray& grids = acoefs->boxArray();
    
    soln.setVal(0.0);

    FArrayBox fab;

    const Real* dx = geom.CellSize();

    const int bho = (max_bndry_order > 2) ? 1 : 0;
    
    for (MFIter mfi(soln); mfi.isValid(); ++mfi)
    {
        int i = mfi.index();
        const Box &reg = grids[i];

        FArrayBox *xfab;
        if (soln.nGrow() == 0) { // need a temporary if soln is the wrong size
            xfab = &soln[mfi];
        } else {
            fab.resize(reg);
            xfab = &fab;
            xfab->copy(soln[mfi], 0, 0, 1);
        }

        HYPRE_StructVectorSetBoxValues(x, Hypre::loV(reg), Hypre::hiV(reg),
                                       xfab->dataPtr());

        FArrayBox *bfab;
        if (rhs.nGrow() == 0) {
            bfab = const_cast<FArrayBox*>(&(rhs[mfi]));
        } else {
            fab.resize(reg);
            bfab->copy(rhs[mfi], 0, 0, 1);
        }
        
        int volume = reg.numPts();
        Real* mat = hypre_CTAlloc(double, stencil_size*volume);

        // build matrix interior
        amrex_hpacoef(BL_TO_FORTRAN_BOX(reg),
                      mat,
                      BL_TO_FORTRAN_ANYD((*acoefs)[mfi]),
                      &scalar_a);

        for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
            amrex_hpbcoef(BL_TO_FORTRAN_BOX(reg),
                          mat,
                          BL_TO_FORTRAN_ANYD((*bcoefs[idim])[mfi]),
                          &scalar_b, dx, &idim);
        }

        const Vector< Vector<BoundCond> > & bcs_i = bndry.bndryConds(i);
        const BndryData::RealTuple      & bcl_i = bndry.bndryLocs(i);

        // add b.c.'s for A matrix and b vector
        for (OrientationIter oit; oit; oit++)
        {
            Orientation ori = oit();
            int cdir(ori);
            int idim = ori.coordDir();
            const int bctype = bcs_i[cdir][0];
            const Real  &bcl = bcl_i[cdir];
            const Mask  &msk = bndry.bndryMasks(ori)[mfi];
      
            amrex_hpmat(BL_TO_FORTRAN_BOX(reg),
                        mat,
                        BL_TO_FORTRAN_ANYD((*bcoefs[idim])[mfi]),
                        BL_TO_FORTRAN_ANYD(msk),
                        &scalar_b, dx, &cdir, &bctype, &bcl, &bho);
        }

        // initialize matrix & vectors
        Array<int,stencil_size> stencil_indices;
        std::iota(stencil_indices.begin(), stencil_indices.end(), 0);
        HYPRE_StructMatrixSetBoxValues(A, Hypre::loV(reg), Hypre::hiV(reg),
                                       stencil_size, stencil_indices.data(), mat);
        HYPRE_StructVectorSetBoxValues(b, Hypre::loV(reg), Hypre::hiV(reg),
                                       bfab->dataPtr());

        hypre_TFree(mat);
    }

    HYPRE_StructMatrixAssemble(A);
    HYPRE_StructVectorAssemble(x); 
    HYPRE_StructVectorAssemble(b); 

    HYPRE_StructPFMGCreate(comm, &solver);
    HYPRE_StructPFMGSetMaxIter(solver, maxiter);
    HYPRE_StructPFMGSetTol(solver, reltol);
    int logging = (verbose >= 2) ? 1 : 0;
    HYPRE_StructPFMGSetLogging(solver, 1);
    HYPRE_StructPFMGSetup(solver, A, b, x);

    if (abstol > 0.0)
    {
        Real bnorm;
        bnorm = hypre_StructInnerProd((hypre_StructVector *) b,
                                      (hypre_StructVector *) b);
        bnorm = std::sqrt(bnorm);

        Real volume = grids.d_numPts();

        Real reltol_new = (bnorm > 0.0
                           ? abstol / bnorm * std::sqrt(volume)
                           : reltol);

        if (reltol_new > reltol) {
            HYPRE_StructPFMGSetTol(solver, reltol_new);
        }
    }

    HYPRE_StructPFMGSolve(solver, A, b, x);

    for (MFIter mfi(soln); mfi.isValid(); ++mfi)
    {
        const Box &reg = mfi.validbox();

        FArrayBox *f;
        if (soln.nGrow() == 0) { // need a temporary if soln is the wrong size
            f = &soln[mfi];
        }
        else {
            fab.resize(reg);
            f = &fab;
        }
        
        Real* vec = f->dataPtr();
        HYPRE_StructVectorGetBoxValues(x, Hypre::loV(reg), Hypre::hiV(reg), vec);
        if (soln.nGrow() != 0) {
            soln[mfi].copy(*f, 0, 0, 1);
        }
    }

    if (verbose >= 2)
    {
        int num_iterations;
        Real res;
        HYPRE_StructPFMGGetNumIterations(solver, &num_iterations);
        HYPRE_StructPFMGGetFinalRelativeResidualNorm(solver, &res);
        
        amrex::Print() << "\n" << num_iterations
                       << " Hypre Multigrid Iterations, Relative Residual "
                       << res << std::endl;
    }

    HYPRE_StructPFMGDestroy(solver);
}

}
