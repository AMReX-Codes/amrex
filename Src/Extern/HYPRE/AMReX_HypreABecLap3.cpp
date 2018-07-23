#include <AMReX_Hypre.H>
#include <AMReX_HypreABec_F.H>
#include <AMReX_VisMF.H>

#ifdef AMREX_USE_EB
#include <AMReX_EBFabFactory.H>
#endif

#include <cmath>
#include <numeric>
#include <limits>
#include <algorithm>
#include <type_traits>

namespace amrex {

constexpr HYPRE_Int HypreABecLap3::regular_stencil_size;
constexpr HYPRE_Int HypreABecLap3::eb_stencil_size;

HypreABecLap3::HypreABecLap3 (const BoxArray& grids,
                              const DistributionMapping& dmap,
                              const Geometry& geom_,
                              MPI_Comm comm_)
    : comm(comm_),
      geom(geom_)
{
    static_assert(AMREX_SPACEDIM > 1, "HypreABecLap2: 1D not supported");

    const int ncomp = 1;
    int ngrow = 0;
    acoefs.define(grids, dmap, ncomp, ngrow);
    acoefs.setVal(0.0);

    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        BoxArray edge_boxes(grids);
        edge_boxes.surroundingNodes(i);
        bcoefs[i].define(edge_boxes, dmap, ncomp, ngrow);
        bcoefs[i].setVal(0.0);
    }
}
    
HypreABecLap3::~HypreABecLap3 ()
{
    HYPRE_IJMatrixDestroy(A);
    A= NULL;
    HYPRE_IJVectorDestroy(b);
    b = NULL;
    HYPRE_IJVectorDestroy(x);
    x = NULL;
    HYPRE_BoomerAMGDestroy(solver);
    solver = NULL;
    m_factory = nullptr;
    m_bndry = nullptr;
    m_maxorder = -1;
}

void
HypreABecLap3::setScalars (Real sa, Real sb)
{
    scalar_a = sa;
    scalar_b = sb;
}

void
HypreABecLap3::setACoeffs (const MultiFab& alpha)
{
    MultiFab::Copy(acoefs, alpha, 0, 0, 1, 0);
}

void
HypreABecLap3::setBCoeffs (const Array<const MultiFab*, AMREX_SPACEDIM>& beta)
{
    for (int idim=0; idim < AMREX_SPACEDIM; idim++) {
        MultiFab::Copy(bcoefs[idim], *beta[idim], 0, 0, 1, 0);
    }
}

void
HypreABecLap3::setVerbose (int _verbose)
{
    verbose = _verbose;
}

void
HypreABecLap3::solve (MultiFab& soln, const MultiFab& rhs,
                      Real rel_tol, Real abs_tol,
                      int max_iter, const BndryData& bndry, int max_bndry_order)
{
    if (solver == NULL || m_bndry != &bndry || m_maxorder != max_bndry_order)
    {
        m_bndry = &bndry;
        m_maxorder = max_bndry_order;
        m_factory = &(rhs.Factory());
        prepareSolver();
    } else {
        HYPRE_IJVectorInitialize(b);
        HYPRE_IJVectorInitialize(x);
    }

#ifdef AMREX_USE_EB
    EBFArrayBoxFactory* ebfactory = dynamic_cast<EBFArrayBoxFactory*>(m_factory);
    const FabArray<EBCellFlagFab>* flags = (ebfactory) ? &(ebfactory->getMultiEBCellFlagFab()) : nullptr;
#endif

    soln.setVal(0.0);
    
    FArrayBox rfab;
    BaseFab<HYPRE_Int> ifab;
    for (MFIter mfi(soln); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();
        const HYPRE_Int nrows = ncells_grid[mfi];
        
#ifdef AMREX_USE_EB
        auto fabtyp = (flags) ? (*flags)[mfi].getType(bx) : FabType::regular;
        if (fabtyp == FabType::covered)
        {
        }
        else if (fabtyp == FabType::singlevalued)
        {
        }
        else
#endif
        {
            ifab.resize(bx);
            ifab.copy(cell_id[mfi],bx);
            
            FArrayBox *xfab;
            if (soln.nGrow() == 0) {
                xfab = &soln[mfi];
            } else {
                xfab = &rfab;
                xfab->resize(bx);
                xfab->setVal(0.0);
            }

            HYPRE_IJVectorSetValues(x, nrows, ifab.dataPtr(), xfab->dataPtr());

            FArrayBox* bfab;
            if (rhs.nGrow() == 0) {
                bfab = const_cast<FArrayBox*>(&rhs[mfi]);
            } else {
                bfab = &rfab;
                bfab->resize(bx);
                bfab->copy(rhs[mfi],bx);
            }

            HYPRE_IJVectorSetValues(b, nrows, ifab.dataPtr(), bfab->dataPtr());
        }
    }

    HYPRE_IJVectorAssemble(x);
    HYPRE_IJVectorAssemble(b);

    HYPRE_ParCSRMatrix par_A = NULL;
    HYPRE_ParVector par_b = NULL;
    HYPRE_ParVector par_x = NULL;
    HYPRE_IJMatrixGetObject(A, (void**)  &par_A);
    HYPRE_IJVectorGetObject(b, (void **) &par_b);
    HYPRE_IJVectorGetObject(x, (void **) &par_x);

    HYPRE_BoomerAMGSetup(solver, par_A, par_b, par_x);

    HYPRE_BoomerAMGSetMinIter(solver, 1);
    HYPRE_BoomerAMGSetMaxIter(solver, max_iter);
    HYPRE_BoomerAMGSetTol(solver, rel_tol);
    if (abs_tol > 0.0)
    {
        Real bnorm = hypre_ParVectorInnerProd(par_b, par_b);
        bnorm = std::sqrt(bnorm);
        
        const BoxArray& grids = rhs.boxArray();
        Real volume = grids.numPts();
        Real rel_tol_new = bnorm > 0.0 ?
            abs_tol / bnorm * std::sqrt(volume) : rel_tol;

        if (rel_tol_new > rel_tol) {
            HYPRE_BoomerAMGSetTol(solver, rel_tol_new);
        }
    }

    HYPRE_BoomerAMGSolve(solver, par_A, par_b, par_x);

    if (verbose >= 2)
    {
        HYPRE_Int num_iterations;
        Real res;
        HYPRE_BoomerAMGGetNumIterations(solver, &num_iterations);
        HYPRE_BoomerAMGGetFinalRelativeResidualNorm(solver, &res);

        amrex::Print() <<"\n" <<  num_iterations
                       << " Hypre BoomerAMG Iterations, Relative Residual "
                       << res << std::endl;
    }

    getSolution(soln);
}

void
HypreABecLap3::getSolution (MultiFab& soln)
{

#ifdef AMREX_USE_EB
    EBFArrayBoxFactory* ebfactory = dynamic_cast<EBFArrayBoxFactory*>(m_factory);
    const FabArray<EBCellFlagFab>* flags = (ebfactory) ? &(ebfactory->getMultiEBCellFlagFab()) : nullptr;
#endif

    FArrayBox rfab;
    BaseFab<HYPRE_Int> ifab;
    for (MFIter mfi(soln); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();
        const HYPRE_Int nrows = ncells_grid[mfi];
        
#ifdef AMREX_USE_EB
        auto fabtyp = (flags) ? (*flags)[mfi].getType(bx) : FabType::regular;
        if (fabtyp == FabType::covered)
        {
        }
        else if (fabtyp == FabType::singlevalued)
        {
        }
        else
#endif
        {
            ifab.resize(bx);
            ifab.copy(cell_id[mfi],bx);
            
            FArrayBox *xfab;
            if (soln.nGrow() == 0) {
                xfab = &soln[mfi];
            } else {
                xfab = &rfab;
                xfab->resize(bx);
            }

            HYPRE_IJVectorGetValues(x, nrows, ifab.dataPtr(), xfab->dataPtr());

            if (soln.nGrow() != 0) {
                soln[mfi].copy(*xfab,bx);
            }
        }
    }
}
   
void
HypreABecLap3::prepareSolver ()
{
    int num_procs, myid;
    MPI_Comm_size(comm, &num_procs);
    MPI_Comm_rank(comm, &myid);

    const BoxArray& ba = acoefs.boxArray();
    const DistributionMapping& dm = acoefs.DistributionMap();
    
#if defined(AMREX_DEBUG) || defined(AMREX_TESTING)
    if (sizeof(HYPRE_Int) < sizeof(long)) {
        long ncells_grids = ba.numPts();
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ncells_grids < static_cast<long>(std::numeric_limits<HYPRE_Int>::max()),
                                         "You might need to configure Hypre with --enable-bigint");
    }
#endif
    
    // how many non-covered cells do we have?
    ncells_grid.define(ba,dm);
    cell_id.define(ba,dm,1,1);

#ifdef AMREX_USE_EB
    EBFArrayBoxFactory* ebfactory = dynamic_cast<EBFArrayBoxFactory*>(m_factory);
    const FabArray<EBCellFlagFab>* flags = (ebfactory) ? &(ebfactory->getMultiEBCellFlagFab()) : nullptr;
#endif

    HYPRE_Int ncells_proc = 0;
#ifdef _OPENMP
#pragma omp parallel reduction(+:ncells_proc)
#endif
    {  BaseFab<HYPRE_Int> ifab;
    for (MFIter mfi(cell_id); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();
        BaseFab<HYPRE_Int>& cid_fab = cell_id[mfi];
        cid_fab.setVal(std::numeric_limits<HYPRE_Int>::lowest());
#ifdef AMREX_USE_EB
        auto fabtyp = (flags) ? (*flags)[mfi].getType(bx) : FabType::regular;
        if (fabtyp == FabType::covered)
        {
            ncells_grid[mfi] = 0;
        }
        else if (fabtyp == FabType::singlevalued)
        {
        }
        else
#endif
        {
            long npts = bx.numPts();
            ncells_grid[mfi] = npts;
            ncells_proc += npts;

            ifab.resize(bx);
            HYPRE_Int* p = ifab.dataPtr();
            for (long i = 0; i < npts; ++i) {
                *p++ = i;
            }
            cid_fab.copy(ifab,bx);
        }
    }}

    Vector<HYPRE_Int> ncells_allprocs(num_procs);
    MPI_Allgather(&ncells_proc, sizeof(HYPRE_Int), MPI_CHAR,
                  ncells_allprocs.data(), sizeof(HYPRE_Int), MPI_CHAR,
                  comm);
    HYPRE_Int proc_begin = 0;
    for (int i = 0; i < myid; ++i) {
        proc_begin += ncells_allprocs[i];
    }

    LayoutData<HYPRE_Int> offset(ba,dm);
    HYPRE_Int proc_end = proc_begin;
    for (MFIter mfi(ncells_grid); mfi.isValid(); ++mfi)
    {
        offset[mfi] = proc_end;
        proc_end += ncells_grid[mfi];
    }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(proc_end == proc_begin+ncells_proc,
                                     "HypreABecLap3::prepareSolver: how did this happend?");

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(cell_id,true); mfi.isValid(); ++mfi)
    {
        cell_id[mfi].plus(offset[mfi], mfi.tilebox());
    }    

    cell_id.FillBoundary(geom.periodicity());

    // Create and initialize A, b & x
    HYPRE_Int ilower = proc_begin;
    HYPRE_Int iupper = proc_end-1;

    //
    HYPRE_IJMatrixCreate(comm, ilower, iupper, ilower, iupper, &A);
    HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);
    HYPRE_IJMatrixInitialize(A);
    //
    HYPRE_IJVectorCreate(comm, ilower, iupper, &b);
    HYPRE_IJVectorSetObjectType(b, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(b);
    //
    HYPRE_IJVectorCreate(comm, ilower, iupper, &x);
    HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(x);
    
    // A.SetValues() & A.assemble()
    const Real* dx = geom.CellSize();
    const int bho = (m_maxorder > 2) ? 1 : 0;
    FArrayBox rfab;
    BaseFab<HYPRE_Int> ifab;
    for (MFIter mfi(acoefs); mfi.isValid(); ++mfi)
    {
        const int i = mfi.index();
        const Box& bx = mfi.validbox();
        
#ifdef AMREX_USE_EB
        const auto fabtyp = (flags) ? (*flags)[mfi].getType(bx) : FabType::regular;
#else
        const auto fabtyp = FabType::regular;
#endif
        if (fabtyp != FabType::covered)
        {
            const HYPRE_Int max_stencil_size = (fabtyp == FabType::regular) ?
                regular_stencil_size : eb_stencil_size;

            ifab.resize(bx,(max_stencil_size+2));
            rfab.resize(bx,max_stencil_size);
            
            const HYPRE_Int nrows = ncells_grid[mfi];
            HYPRE_Int* ncols = ifab.dataPtr(0);
            HYPRE_Int* rows  = ifab.dataPtr(1);
            HYPRE_Int* cols  = ifab.dataPtr(2);
            Real*      mat   = rfab.dataPtr();

            Array<int,AMREX_SPACEDIM*2> bctype;
            Array<Real,AMREX_SPACEDIM*2> bcl;
            const Vector< Vector<BoundCond> > & bcs_i = m_bndry->bndryConds(i);
            const BndryData::RealTuple        & bcl_i = m_bndry->bndryLocs(i);
            for (OrientationIter oit; oit; oit++) {
                int cdir(oit());
                bctype[cdir] = bcs_i[cdir][0];
                bcl[cdir]  = bcl_i[cdir];
            }
            
            if (fabtyp == FabType::regular)
            {
                amrex_hpijmatrix(BL_TO_FORTRAN_BOX(bx),
                                 &nrows, ncols, rows, cols, mat,
                                 BL_TO_FORTRAN_ANYD(cell_id[mfi]),
                                 &(offset[mfi]),
                                 BL_TO_FORTRAN_ANYD(acoefs[mfi]),
                                 AMREX_D_DECL(BL_TO_FORTRAN_ANYD(bcoefs[0][mfi]),
                                              BL_TO_FORTRAN_ANYD(bcoefs[1][mfi]),
                                              BL_TO_FORTRAN_ANYD(bcoefs[2][mfi])),
                                 &scalar_a, &scalar_b, dx,
                                 bctype.data(), bcl.data(), &bho);
            }
#ifdef AMREX_USE_EB
            else
            {
                amrex::Print() << "HypreABecLap3: set mat todo\n";
            }
#endif
            HYPRE_IJMatrixSetValues(A,nrows,ncols,rows,cols,mat);
        }
    }
    HYPRE_IJMatrixAssemble(A);

    // Create solver
    HYPRE_BoomerAMGCreate(&solver);

#if 0
    // Set some parameters (See Reference Manual for more parameters)
    // Falgout coarsening with modified classical interpolation
    HYPRE_BoomerAMGSetOldDefault(solver);
    HYPRE_BoomerAMGSetCoarsenType(solver, 6);
    HYPRE_BoomerAMGSetCycleType(solver, 1);
    HYPRE_BoomerAMGSetRelaxType(solver, 6);   /* G-S/Jacobi hybrid relaxation */
    HYPRE_BoomerAMGSetRelaxOrder(solver, 1);   /* uses C/F relaxation */
    HYPRE_BoomerAMGSetNumSweeps(solver, 2);   /* Sweeeps on each level */
    HYPRE_BoomerAMGSetMaxLevels(solver, 20);  /* maximum number of levels */
    HYPRE_BoomerAMGSetStrongThreshold(solver, 0.6);
#endif

    int logging = (verbose >= 2) ? 1 : 0;
    HYPRE_BoomerAMGSetLogging(solver, logging);
}
    
}  // namespace amrex
