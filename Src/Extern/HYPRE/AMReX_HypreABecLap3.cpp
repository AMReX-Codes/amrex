#include <AMReX_HypreABecLap3.H>
#include <AMReX_VisMF.H>

#include <AMReX_Habec_K.H>

#ifdef AMREX_USE_EB
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_MultiCutFab.H>
#include <AMReX_EBFabFactory.H>
#endif

#include <cmath>
#include <numeric>
#include <limits>
#include <algorithm>
#include <type_traits>

namespace amrex {

HypreABecLap3::HypreABecLap3 (const BoxArray& grids, const DistributionMapping& dmap,
                              const Geometry& geom_, MPI_Comm comm_)
    : Hypre(grids, dmap, geom_, comm_)
{
}
    
HypreABecLap3::~HypreABecLap3 ()
{
    HYPRE_IJMatrixDestroy(A);
    A = NULL;
    HYPRE_IJVectorDestroy(b);
    b = NULL;
    HYPRE_IJVectorDestroy(x);
    x = NULL;
    HYPRE_BoomerAMGDestroy(solver);
    solver = NULL;
}

void
HypreABecLap3::solve (MultiFab& soln, const MultiFab& rhs, Real rel_tol, Real abs_tol,
                      int max_iter, const BndryData& bndry, int max_bndry_order)
{
    Gpu::LaunchSafeGuard lsg(false); // xxxxx TODO: gpu

    BL_PROFILE("HypreABecLap3::solve()");

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
    
    HYPRE_IJVectorInitialize(b);
    HYPRE_IJVectorInitialize(x);
    //
    loadVectors(soln, rhs);
    //
    HYPRE_IJVectorAssemble(x);
    HYPRE_IJVectorAssemble(b);
    
    HYPRE_ParCSRMatrix par_A = NULL;
    HYPRE_ParVector par_b = NULL;
    HYPRE_ParVector par_x = NULL;
    HYPRE_IJMatrixGetObject(A, (void**)  &par_A);
    HYPRE_IJVectorGetObject(b, (void **) &par_b);
    HYPRE_IJVectorGetObject(x, (void **) &par_x);

    HYPRE_BoomerAMGSetMinIter(solver, 1);
    HYPRE_BoomerAMGSetMaxIter(solver, max_iter);
    HYPRE_BoomerAMGSetTol(solver, rel_tol);
    if (abs_tol > 0.0)
    {
        Real bnorm = hypre_ParVectorInnerProd(par_b, par_b);
        bnorm = std::sqrt(bnorm);
        
        const BoxArray& grids = rhs.boxArray();
        Real volume = grids.numPts();
        Real rel_tol_new = (bnorm > 0.0) ? (abs_tol / bnorm * std::sqrt(volume)) : rel_tol;

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
                       << " Hypre IJ BoomerAMG Iterations, Relative Residual "
                       << res << std::endl;
    }

    getSolution(soln);
}

void
HypreABecLap3::getSolution (MultiFab& soln)
{
#ifdef AMREX_USE_EB
    auto ebfactory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory);
    const FabArray<EBCellFlagFab>* flags = (ebfactory) ? &(ebfactory->getMultiEBCellFlagFab()) : nullptr;
#endif

    FArrayBox rfab;
    for (MFIter mfi(soln); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();
        const HYPRE_Int nrows = ncells_grid[mfi];
        
#ifdef AMREX_USE_EB
        auto fabtyp = (flags) ? (*flags)[mfi].getType(bx) : FabType::regular;
#else
        auto fabtyp = FabType::regular;
#endif
        if (fabtyp != FabType::covered)
        {
            FArrayBox *xfab;
#ifdef AMREX_USE_EB
            if (fabtyp != FabType::regular)
            {
                xfab = &rfab;
                xfab->resize(bx);
            }
            else
#endif
            {
                if (soln.nGrow() == 0) {
                    xfab = &soln[mfi];
                } else {
                    xfab = &rfab;
                    xfab->resize(bx);
                }
            }
                
            HYPRE_IJVectorGetValues(x, nrows, cell_id_vec[mfi].data(), xfab->dataPtr());

            if (fabtyp == FabType::regular && soln.nGrow() != 0)
            {
                soln[mfi].copy<RunOn::Host>(*xfab,bx);
            }
#ifdef AMREX_USE_EB
            else if (fabtyp != FabType::regular)
            {
                amrex_hpeb_copy_from_vec(bx,
                                         soln[mfi],
                                         xfab->dataPtr(),
                                         (*flags)[mfi]);
            }
#endif
        }
    }
}
   
void
HypreABecLap3::prepareSolver ()
{
    BL_PROFILE("HypreABecLap3::prepareSolver()");
    
    int num_procs, myid;
    MPI_Comm_size(comm, &num_procs);
    MPI_Comm_rank(comm, &myid);

    const BoxArray& ba = acoefs.boxArray();
    const DistributionMapping& dm = acoefs.DistributionMap();
    
#if defined(AMREX_DEBUG) || defined(AMREX_TESTING)
    if (sizeof(HYPRE_Int) < sizeof(Long)) {
        Long ncells_grids = ba.numPts();
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ncells_grids < static_cast<Long>(std::numeric_limits<HYPRE_Int>::max()),
                                         "You might need to configure Hypre with --enable-bigint");
    }
#endif

    // how many non-covered cells do we have?
    ncells_grid.define(ba,dm);
    cell_id.define(ba,dm,1,1);
    cell_id_vec.define(ba,dm);

#ifdef AMREX_USE_EB
    auto ebfactory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory);
    const FabArray<EBCellFlagFab>* flags = (ebfactory) ? &(ebfactory->getMultiEBCellFlagFab()) : nullptr;
    const MultiFab* vfrac = (ebfactory) ? &(ebfactory->getVolFrac()) : nullptr;
    auto area = (ebfactory) ? ebfactory->getAreaFrac()
        : Array<const MultiCutFab*,AMREX_SPACEDIM>{AMREX_D_DECL(nullptr,nullptr,nullptr)};
    auto fcent = (ebfactory) ? ebfactory->getFaceCent()
        : Array<const MultiCutFab*,AMREX_SPACEDIM>{AMREX_D_DECL(nullptr,nullptr,nullptr)};
    auto barea = (ebfactory) ? &(ebfactory->getBndryArea()) : nullptr;
    auto bcent = (ebfactory) ? &(ebfactory->getBndryCent()) : nullptr;
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
        cid_fab.setVal<RunOn::Host>(std::numeric_limits<HYPRE_Int>::lowest());
#ifdef AMREX_USE_EB
        auto fabtyp = (flags) ? (*flags)[mfi].getType(bx) : FabType::regular;
        if (fabtyp == FabType::covered)
        {
            ncells_grid[mfi] = 0;
        }
        else if (fabtyp == FabType::singlevalued)
        {
            amrex_hpeb_fill_cellid(bx,
                                   ncells_grid[mfi],
                                   cid_fab,
                                   (*flags)[mfi]);

            ncells_proc += ncells_grid[mfi];
        }
        else
#endif
        {
            Long npts = bx.numPts();
            ncells_grid[mfi] = npts;
            ncells_proc += npts;

            ifab.resize(bx);
            HYPRE_Int* p = ifab.dataPtr();
            for (Long i = 0; i < npts; ++i) {
                *p++ = i;
            }
            cid_fab.copy<RunOn::Host>(ifab,bx);
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

#ifdef AMREX_DEBUG
    HYPRE_Int ncells_total = 0;
    for (auto n : ncells_allprocs) {
        ncells_total += n;
    }
#endif

    LayoutData<HYPRE_Int> offset(ba,dm);
    HYPRE_Int proc_end = proc_begin;
    for (MFIter mfi(ncells_grid); mfi.isValid(); ++mfi)
    {
        offset[mfi] = proc_end;
        proc_end += ncells_grid[mfi];
    }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(proc_end == proc_begin+ncells_proc,
                                     "HypreABecLap3::prepareSolver: how did this happen?");

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(cell_id,true); mfi.isValid(); ++mfi)
    {
        cell_id[mfi].plus<RunOn::Host>(offset[mfi], mfi.tilebox());
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
    //
    HYPRE_IJVectorCreate(comm, ilower, iupper, &x);
    HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR);
    
    // A.SetValues() & A.assemble()

    const Real* dx = geom.CellSize();
    const int bho = (m_maxorder > 2) ? 1 : 0;
    FArrayBox foo(Box::TheUnitBox());
    const int is_eb_dirichlet = m_eb_b_coeffs != nullptr;

    for (MFIter mfi(acoefs); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();
        
#ifdef AMREX_USE_EB
        const auto fabtyp = (flags) ? (*flags)[mfi].getType(bx) : FabType::regular;
#else
        const auto fabtyp = FabType::regular;
#endif
        if (fabtyp != FabType::covered)
        {

            const HYPRE_Int nrows = ncells_grid[mfi];
            cell_id_vec[mfi].resize(nrows);

            Gpu::ManagedDeviceVector<HYPRE_Int> ncolsg(nrows,0);
            Gpu::ManagedDeviceVector<HYPRE_Int> colsg(nrows*(AMREX_SPACEDIM*2+1),0);
            Gpu::ManagedDeviceVector<Real> matg(nrows*(AMREX_SPACEDIM*2+1),0.0);

            GpuArray<int,AMREX_SPACEDIM*2> bctype;
            GpuArray<Real,AMREX_SPACEDIM*2> bcl;
            const Vector< Vector<BoundCond> > & bcs_i = m_bndry->bndryConds(mfi);
            const BndryData::RealTuple        & bcl_i = m_bndry->bndryLocs(mfi);
            for (OrientationIter oit; oit; oit++) {
                int cdir(oit());
                bctype[cdir] = bcs_i[cdir][0];
                bcl[cdir]  = bcl_i[cdir];
            }
            
            if (fabtyp == FabType::regular)
            {
                amrex_hpijmatrix(bx,
                                 nrows, ncolsg.dataPtr(), 
                                 cell_id_vec[mfi].dataPtr(), 
                                 colsg.dataPtr(), matg.dataPtr(),
                                 cell_id[mfi], 
                                 offset[mfi],
                                 diaginv[mfi],
                                 acoefs[mfi],
                                 bcoefs[0][mfi],
                                 bcoefs[1][mfi],
#if (AMREX_SPACEDIM == 3)
                                 bcoefs[2][mfi],
#endif
                                 scalar_a, scalar_b, dx,
                                 bctype, bcl, bho);

            }
#ifdef AMREX_USE_EB
            else
            {
                FArrayBox const& beb = (is_eb_dirichlet) ? (*m_eb_b_coeffs)[mfi] : foo;
                int size_vec = std::pow(3,AMREX_SPACEDIM);
                colsg.resize(nrows*size_vec);
                matg.resize(nrows*size_vec);
                amrex_hpeb_ijmatrix(bx,
                                    nrows, ncolsg.dataPtr(),
                                    cell_id_vec[mfi].dataPtr(),
                                    colsg.dataPtr(), matg.dataPtr(),
                                    cell_id[mfi],
                                    offset[mfi], diaginv[mfi],
                                    acoefs[mfi], bcoefs[0][mfi],
                                    bcoefs[1][mfi], 
#if (AMREX_SPACEDIM == 3)
                                    bcoefs[2][mfi],
#endif
                                    (*flags)[mfi],
                                    (*vfrac)[mfi],
                                    (*area[0])[mfi], (*area[1])[mfi],
#if (AMREX_SPACEDIM == 3)
                                    (*area[2])[mfi],
#endif 
                                    (*fcent[0])[mfi], (*fcent[1])[mfi], 
#if (AMREX_SPACEDIM == 3)
                                    (*fcent[2])[mfi],
#endif
                                    (*barea)[mfi],
                                    (*bcent)[mfi],
                                    beb, is_eb_dirichlet,
                                    scalar_a, scalar_b, dx,
                                    bctype, bcl, bho);

            }
#endif

            HYPRE_Int* rows = cell_id_vec[mfi].data();
            HYPRE_Int* ncols = (HYPRE_Int*) ncolsg.dataPtr();
            HYPRE_Int* cols = (HYPRE_Int*) colsg.dataPtr();
            Real* mat = (Real*) matg.dataPtr();

#ifdef AMREX_DEBUG
            HYPRE_Int nvalues = 0;
            for (HYPRE_Int i = 0; i < nrows; ++i) {
                nvalues += ncols[i];
            }
            for (HYPRE_Int i = 0; i < nvalues; ++i) {
                AMREX_ASSERT(cols[i] >= 0 && cols[i] < ncells_total);
            }
#endif

            HYPRE_IJMatrixSetValues(A,nrows,ncols,rows,cols,mat);
        }
    }
    HYPRE_IJMatrixAssemble(A);

    // Create solver
    HYPRE_BoomerAMGCreate(&solver);

    if (old_default) HYPRE_BoomerAMGSetOldDefault(solver); // Falgout coarsening with modified classical interpolation
//    HYPRE_BoomerAMGSetCoarsenType(solver, 6);
//    HYPRE_BoomerAMGSetCycleType(solver, 1);
    HYPRE_BoomerAMGSetRelaxType(solver, relax_type);
    HYPRE_BoomerAMGSetRelaxOrder(solver, relax_order);
    HYPRE_BoomerAMGSetNumSweeps(solver, num_sweeps);
    HYPRE_BoomerAMGSetStrongThreshold(solver, strong_threshold);

    int logging = (verbose >= 2) ? 1 : 0;
    HYPRE_BoomerAMGSetLogging(solver, logging);

    HYPRE_ParCSRMatrix par_A = NULL;
    HYPRE_IJMatrixGetObject(A, (void**)  &par_A);
    HYPRE_BoomerAMGSetup(solver, par_A, NULL, NULL);
}

void
HypreABecLap3::loadVectors (MultiFab& soln, const MultiFab& rhs)
{
    BL_PROFILE("HypreABecLap3::loadVectors()");
    
#ifdef AMREX_USE_EB
    auto ebfactory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory);
    const FabArray<EBCellFlagFab>* flags = (ebfactory) ? &(ebfactory->getMultiEBCellFlagFab()) : nullptr;
#endif

    soln.setVal(0.0);
    
    FArrayBox vecfab, rhsfab;
    for (MFIter mfi(soln); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();
        const HYPRE_Int nrows = ncells_grid[mfi];
        
#ifdef AMREX_USE_EB
        const auto fabtyp = (flags) ? (*flags)[mfi].getType(bx) : FabType::regular;
#else
        const auto fabtyp = FabType::regular;
#endif
        if (fabtyp != FabType::covered)
        {
            // soln has been set to zero.
            HYPRE_IJVectorSetValues(x, nrows, cell_id_vec[mfi].data(), soln[mfi].dataPtr());

            rhsfab.resize(bx);
            rhsfab.copy<RunOn::Host>(rhs[mfi],bx);
            rhsfab.mult<RunOn::Host>(diaginv[mfi]);
            
            FArrayBox* bfab;                            
#ifdef AMREX_USE_EB
            if (fabtyp != FabType::regular)
            {
                bfab = &vecfab;
                bfab->resize(bx);
                amrex_hpeb_copy_to_vec(bx,
                                       rhsfab,
                                       bfab->dataPtr(),
                                       (*flags)[mfi]);
            }
            else
#endif
            {
                bfab = &rhsfab;
            }

            HYPRE_IJVectorSetValues(b, nrows, cell_id_vec[mfi].data(), bfab->dataPtr());
        }
    }
}
    
}  // namespace amrex
