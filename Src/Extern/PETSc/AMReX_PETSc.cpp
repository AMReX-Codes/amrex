
#include <petscksp.h>
#include <AMReX_PETSc.H>

#include <cmath>
#include <numeric>
#include <limits>
#include <algorithm>
#include <type_traits>

namespace amrex {

constexpr PetscInt PETScABecLap::regular_stencil_size;
constexpr PetscInt PETScABecLap::eb_stencil_size;

PETScABecLap::PETScABecLap (const BoxArray& grids, const DistributionMapping& dmap,
                            const Geometry& geom_, MPI_Comm comm_)
    : comm(comm_),
      geom(geom_)
{
    static_assert(AMREX_SPACEDIM > 1, "PETScABecLap: 1D not supported");
    static_assert(std::is_same<Real, PetscScalar>::value, "amrex::Real != PetscScalar");
    static_assert(std::is_same<HYPRE_Int, PetscInt>::value, "HYPRE_Int != PetscInt");

    const int ncomp = 1;
    int ngrow = 0;
    acoefs.define(grids, dmap, ncomp, ngrow);
    acoefs.setVal(0.0);
    
#ifdef AMREX_USE_EB
    ngrow = 1;
#endif
    
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        BoxArray edge_boxes(grids);
        edge_boxes.surroundingNodes(i);
        bcoefs[i].define(edge_boxes, dmap, ncomp, ngrow);
        bcoefs[i].setVal(0.0);
    }

    diaginv.define(grids,dmap,ncomp,0);

    PETSC_COMM_WORLD = comm;
    PetscInitialize(0, 0, 0, 0);
}

PETScABecLap::~PETScABecLap ()
{
    MatDestroy(&A);
    A = nullptr;

    VecDestroy(&b);
    b = nullptr;

    VecDestroy(&x);
    x = nullptr;

    KSPDestroy(&solver);
    solver = nullptr;

    m_factory = nullptr;
    m_bndry = nullptr;
    m_maxorder = -1;
    
    PetscFinalize();
}

void
PETScABecLap::setScalars (Real sa, Real sb)
{
    scalar_a = sa;
    scalar_b = sb;
}

void
PETScABecLap::setACoeffs (const MultiFab& alpha)
{
    MultiFab::Copy(acoefs, alpha, 0, 0, 1, 0);
}

void
PETScABecLap::setBCoeffs (const Array<const MultiFab*, BL_SPACEDIM>& beta)
{
    for (int idim=0; idim < AMREX_SPACEDIM; idim++) {
        const int ng = std::min(bcoefs[idim].nGrow(), beta[idim]->nGrow());
        MultiFab::Copy(bcoefs[idim], *beta[idim], 0, 0, 1, ng);
    }
}

void
PETScABecLap::setVerbose (int _verbose)
{
    verbose = _verbose;
}

void
PETScABecLap::solve (MultiFab& soln, const MultiFab& rhs, Real rel_tol, Real abs_tol, 
                     int max_iter, const BndryData& bndry, int max_bndry_order)
{
    BL_PROFILE("PETScABecLap::solve()");

    if (solver == nullptr || m_bndry != &bndry || m_maxorder != max_bndry_order)
    {
        m_bndry = &bndry;
        m_maxorder = max_bndry_order;
        m_factory = &(rhs.Factory());
        prepareSolver();
    }

    loadVectors(soln, rhs);

    KSPSetTolerances(solver, rel_tol, PETSC_DEFAULT, PETSC_DEFAULT, max_iter);
    
    KSPSolve(solver, b, x);

    if (verbose >= 2)
    {
        PetscInt niters;
        Real res;
        KSPGetIterationNumber(solver, &niters);
        KSPGetResidualNorm(solver, &res);
        amrex::Print() <<"\n" <<  niters << " PETSc Iterations, Residual Norm " << res << std::endl;
    }

    getSolution(soln);
}

void
PETScABecLap::prepareSolver ()
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
                                         "PetscInt is too short");
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
            amrex_hpeb_fill_cellid(BL_TO_FORTRAN_BOX(bx),
                                   &ncells_grid[mfi],
                                   BL_TO_FORTRAN_ANYD(cid_fab),
                                   BL_TO_FORTRAN_ANYD((*flags)[mfi]));
            ncells_proc += ncells_grid[mfi];
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


    KSPCreate(comm, &solver);

    // create A

    KSPSetOperators(solver, A, A);
    PC pc;
    KSPGetPC(solver, &pc);
    PCSetType(pc, PCGAMG);
    KSPSetUp(solver);
}

void
PETScABecLap::loadVectors (MultiFab& soln, const MultiFab& rhs)
{
    soln.setVal(0.0);
}

void
PETScABecLap::getSolution (MultiFab& soln)
{
}

}
