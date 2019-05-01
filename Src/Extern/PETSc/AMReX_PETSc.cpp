
#include <petscksp.h>
#include <AMReX_PETSc.H>

#ifdef AMREX_USE_EB
#include <AMReX_MultiCutFab.H>
#include <AMReX_EBFabFactory.H>
#endif

#include <AMReX_HypreABec_F.H>
#include <cmath>
#include <numeric>
#include <limits>
#include <algorithm>
#include <type_traits>

#ifdef AMREX_USE_HYPRE
#include "_hypre_utilities.h"
#endif


namespace amrex {

constexpr PetscInt PETScABecLap::regular_stencil_size;
constexpr PetscInt PETScABecLap::eb_stencil_size;


std::unique_ptr<PETScABecLap>
makePetsc (const BoxArray& grids, const DistributionMapping& dmap,
           const Geometry& geom, MPI_Comm comm_)
{
        return std::unique_ptr<PETScABecLap>(new PETScABecLap(grids, dmap, geom, comm_));
}


PETScABecLap::PETScABecLap (const BoxArray& grids, const DistributionMapping& dmap,
                            const Geometry& geom_, MPI_Comm comm_)
    : comm(comm_),
      geom(geom_)
{
    static_assert(AMREX_SPACEDIM > 1, "PETScABecLap: 1D not supported");
    static_assert(std::is_same<Real, PetscScalar>::value, "amrex::Real != PetscScalar");
#ifdef AMREX_USE_HYPRE
    static_assert(std::is_same<HYPRE_Int, PetscInt>::value, "HYPRE_Int != PetscInt");
#endif

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
    else
    {
        m_factory = &(rhs.Factory());
    }

    loadVectors(soln, rhs);
    //
    VecAssemblyBegin(x); 
    VecAssemblyEnd(x); 
    //
    VecAssemblyBegin(b); 
    VecAssemblyEnd(b); 
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
    MPI_Comm_size(PETSC_COMM_WORLD, &num_procs);
    MPI_Comm_rank(PETSC_COMM_WORLD, &myid);

    const BoxArray& ba = acoefs.boxArray();
    const DistributionMapping& dm = acoefs.DistributionMap();
    
#if defined(AMREX_DEBUG) || defined(AMREX_TESTING)
    if (sizeof(PetscInt) < sizeof(long)) {
        long ncells_grids = ba.numPts();
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ncells_grids < static_cast<long>(std::numeric_limits<PetscInt>::max()),
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
    auto barea = (ebfactory) ? &(ebfactory->getBndryArea()) : nullptr;
    auto bcent = (ebfactory) ? &(ebfactory->getBndryCent()) : nullptr;
#endif

    PetscInt ncells_proc = 0;
#ifdef _OPENMP
#pragma omp parallel reduction(+:ncells_proc)
#endif
    {  BaseFab<PetscInt> ifab;
    for (MFIter mfi(cell_id); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();
        BaseFab<PetscInt>& cid_fab = cell_id[mfi];
        cid_fab.setVal(std::numeric_limits<PetscInt>::lowest());
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
            PetscInt* p = ifab.dataPtr();
            for (long i = 0; i < npts; ++i) {
                *p++ = i;
            }
            cid_fab.copy(ifab,bx);
        }
    }
    }

    Vector<PetscInt> ncells_allprocs(num_procs);
    MPI_Allgather(&ncells_proc, sizeof(PetscInt), MPI_CHAR,
                  ncells_allprocs.data(), sizeof(PetscInt), MPI_CHAR,
                  PETSC_COMM_WORLD);
    PetscInt proc_begin = 0;
    for (int i = 0; i < myid; ++i) {
        proc_begin += ncells_allprocs[i];
    }
    PetscInt ncells_world = 0;
    for (auto i : ncells_allprocs) {
        ncells_world += i;
    }

    LayoutData<PetscInt> offset(ba,dm);
    PetscInt proc_end = proc_begin;
    for (MFIter mfi(ncells_grid); mfi.isValid(); ++mfi)
    {
        offset[mfi] = proc_end;
        proc_end += ncells_grid[mfi];
    }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(proc_end == proc_begin+ncells_proc,
                                     "PETScABecLap::prepareSolver: how did this happen?");

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(cell_id,true); mfi.isValid(); ++mfi)
    {
        cell_id[mfi].plus(offset[mfi], mfi.tilebox());
    }    

    cell_id.FillBoundary(geom.periodicity());

    // estimated amount of block diag elements
    PetscInt d_nz = (eb_stencil_size + regular_stencil_size) / 2;
    // estimated amount of block off diag elements
    PetscInt o_nz  = d_nz / 2;
    MatCreate(PETSC_COMM_WORLD, &A); 
    MatSetType(A, MATMPIAIJ);
    MatSetSizes(A, ncells_proc, ncells_proc, ncells_world, ncells_world); 
    MatMPIAIJSetPreallocation(A, d_nz, NULL, o_nz, NULL );
    //Maybe an over estimate of the diag/off diag #of non-zero entries, so we turn off malloc warnings
    MatSetUp(A); 
    MatSetOption(A, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_FALSE); 
 
    // A.SetValues
    const Real* dx = geom.CellSize();
    const int bho = (m_maxorder > 2) ? 1 : 0;
    FArrayBox rfab;
    BaseFab<PetscInt> ifab;
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
            const PetscInt max_stencil_size = (fabtyp == FabType::regular) ?
                regular_stencil_size : eb_stencil_size;

            ifab.resize(bx,(max_stencil_size+1));
            rfab.resize(bx,max_stencil_size);

            const PetscInt nrows = ncells_grid[mfi];
            cell_id_vec[mfi].resize(nrows);
            PetscInt* rows = cell_id_vec[mfi].data();
            PetscInt* ncols = ifab.dataPtr(0);
            PetscInt* cols  = ifab.dataPtr(1);
            Real*      mat   = rfab.dataPtr();

            Array<int,AMREX_SPACEDIM*2> bctype;
            Array<Real,AMREX_SPACEDIM*2> bcl;
            const Vector< Vector<BoundCond> > & bcs_i = m_bndry->bndryConds(mfi);
            const BndryData::RealTuple        & bcl_i = m_bndry->bndryLocs(mfi);
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
                                 BL_TO_FORTRAN_ANYD(diaginv[mfi]),
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
                FArrayBox const& beb = (is_eb_dirichlet) ? (*m_eb_b_coeffs)[mfi] : foo;
                
                amrex_hpeb_ijmatrix(BL_TO_FORTRAN_BOX(bx),
                                    &nrows, ncols, rows, cols, mat,
                                    BL_TO_FORTRAN_ANYD(cell_id[mfi]),
                                    &(offset[mfi]),
                                    BL_TO_FORTRAN_ANYD(diaginv[mfi]),
                                    BL_TO_FORTRAN_ANYD(acoefs[mfi]),
                                    AMREX_D_DECL(BL_TO_FORTRAN_ANYD(bcoefs[0][mfi]),
                                                 BL_TO_FORTRAN_ANYD(bcoefs[1][mfi]),
                                                 BL_TO_FORTRAN_ANYD(bcoefs[2][mfi])),
                                    BL_TO_FORTRAN_ANYD((*flags)[mfi]),
                                    BL_TO_FORTRAN_ANYD((*vfrac)[mfi]),
                                    AMREX_D_DECL(BL_TO_FORTRAN_ANYD((*area[0])[mfi]),
                                                 BL_TO_FORTRAN_ANYD((*area[1])[mfi]),
                                                 BL_TO_FORTRAN_ANYD((*area[2])[mfi])),
                                    AMREX_D_DECL(BL_TO_FORTRAN_ANYD((*fcent[0])[mfi]),
                                                 BL_TO_FORTRAN_ANYD((*fcent[1])[mfi]),
                                                 BL_TO_FORTRAN_ANYD((*fcent[2])[mfi])),
                                    BL_TO_FORTRAN_ANYD((*barea)[mfi]),
                                    BL_TO_FORTRAN_ANYD((*bcent)[mfi]),
                                    BL_TO_FORTRAN_ANYD(beb), &is_eb_dirichlet,
                                    &scalar_a, &scalar_b, dx,
                                    bctype.data(), bcl.data(), &bho);
            }
#endif
            //Load in by row! 
            int matid = 0; 
            for (int rit = 0; rit < nrows; ++rit)
            {
                MatSetValues(A, 1, &rows[rit], ncols[rit], &cols[matid], &mat[matid], INSERT_VALUES);
                matid += ncols[rit];
            }
        }
    }

    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    // create solver
    KSPCreate(PETSC_COMM_WORLD, &solver);
    KSPSetOperators(solver, A, A);

    // Set up preconditioner
    PC pc;
    KSPGetPC(solver, &pc);

    // Classic AMG
    PCSetType(pc, PCGAMG);
    PCGAMGSetType(pc, PCGAMGAGG);
    PCGAMGSetNSmooths(pc,0); 
//    PCSetType(pc, PCJACOBI); 

    
// we are not using command line options    KSPSetFromOptions(solver);
    // create b & x
    VecCreateMPI(PETSC_COMM_WORLD, ncells_proc, ncells_world, &x);
    VecDuplicate(x, &b);
}

void
PETScABecLap::loadVectors (MultiFab& soln, const MultiFab& rhs)
{
    BL_PROFILE("PETScABecLap::loadVectors()");

#ifdef AMREX_USE_EB
    auto ebfactory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory);
    const FabArray<EBCellFlagFab>* flags = (ebfactory) ? &(ebfactory->getMultiEBCellFlagFab()) : nullptr;
#endif
    soln.setVal(0.0);

    FArrayBox vecfab, rhsfab;
    for (MFIter mfi(soln); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();
        const PetscInt nrows = ncells_grid[mfi];

#ifdef AMREX_USE_EB
        const auto fabtyp = (flags) ? (*flags)[mfi].getType(bx) : FabType::regular;
#else
        const auto fabtyp = FabType::regular;
#endif
        if (fabtyp != FabType::covered)
        {
            // soln has been set to zero.
            VecSetValues(x, nrows, cell_id_vec[mfi].data(), soln[mfi].dataPtr(), INSERT_VALUES); 
            rhsfab.resize(bx);
            rhsfab.copy(rhs[mfi],bx);
            rhsfab.mult(diaginv[mfi]);

            FArrayBox* bfab;
#ifdef AMREX_USE_EB
            if (fabtyp != FabType::regular)
            {
                bfab = &vecfab;
                bfab->resize(bx);
                amrex_hpeb_copy_to_vec(BL_TO_FORTRAN_BOX(bx),
                                       BL_TO_FORTRAN_ANYD(rhsfab),
                                       bfab->dataPtr(), &nrows,
                                       BL_TO_FORTRAN_ANYD((*flags)[mfi])); // */
            }
            else
#endif
            {
                bfab = &rhsfab;
            }
            VecSetValues(b, nrows, cell_id_vec[mfi].data(), bfab->dataPtr(), INSERT_VALUES); 
        }
    }
}

void
PETScABecLap::getSolution (MultiFab& soln)
{
#ifdef AMREX_USE_EB
    auto ebfactory = dynamic_cast<EBFArrayBoxFactory const*>(m_factory);
    const FabArray<EBCellFlagFab>* flags = (ebfactory) ? &(ebfactory->getMultiEBCellFlagFab()) : nullptr;
#endif

    FArrayBox rfab;
    for (MFIter mfi(soln); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();
        const PetscInt nrows = ncells_grid[mfi];

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

            VecGetValues(x, nrows, cell_id_vec[mfi].data(), xfab->dataPtr()); 
            if (fabtyp == FabType::regular && soln.nGrow() != 0)
            {
                soln[mfi].copy(*xfab,bx);
            }
#ifdef AMREX_USE_EB
            else if (fabtyp != FabType::regular)
            {
                amrex_hpeb_copy_from_vec(BL_TO_FORTRAN_BOX(bx),
                                         BL_TO_FORTRAN_ANYD(soln[mfi]),
                                         xfab->dataPtr(), &nrows,
                                         BL_TO_FORTRAN_ANYD((*flags)[mfi])); // */
            }
#endif
        }
    }
}

}
