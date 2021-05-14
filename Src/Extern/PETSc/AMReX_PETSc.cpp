
#include <petscksp.h>
#include <AMReX_PETSc.H>

#ifdef AMREX_USE_EB
#include <AMReX_MultiCutFab.H>
#include <AMReX_EBFabFactory.H>
#endif

#include <AMReX_Habec_K.H>

#include <cmath>
#include <numeric>
#include <limits>
#include <algorithm>
#include <type_traits>

namespace amrex {

constexpr PetscInt PETScABecLap::regular_stencil_size;
constexpr PetscInt PETScABecLap::eb_stencil_size;


std::unique_ptr<PETScABecLap>
makePetsc (const BoxArray& grids, const DistributionMapping& dmap,
           const Geometry& geom, MPI_Comm comm_)
{
    return std::make_unique<PETScABecLap>(grids, dmap, geom, comm_);
}


PETScABecLap::PETScABecLap (const BoxArray& grids, const DistributionMapping& dmap,
                            const Geometry& geom_, MPI_Comm comm_)
    : comm(comm_),
      geom(geom_)
{
    static_assert(AMREX_SPACEDIM > 1, "PETScABecLap: 1D not supported");
    static_assert(std::is_same<Real, PetscScalar>::value, "amrex::Real != PetscScalar");

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
PETScABecLap::solve (MultiFab& soln, const MultiFab& rhs, Real rel_tol, Real /*abs_tol*/,
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
    if (sizeof(PetscInt) < sizeof(Long)) {
        Long ncells_grids = ba.numPts();
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ncells_grids < static_cast<Long>(std::numeric_limits<PetscInt>::max()),
                                         "PetscInt is too short");
    }
#endif

    static_assert(std::is_signed<PetscInt>::value, "PetscInt is assumed to be signed");

    // how many non-covered cells do we have?
    ncells_grid.define(ba,dm);
    cell_id.define(ba,dm,1,1);
    cell_id_vec.define(ba,dm,1,0);

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
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion()) reduction(+:ncells_proc)
#endif
    for (MFIter mfi(cell_id); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();
        const Box& gbx = amrex::grow(bx,1);
        Array4<PetscInt> const& cid_arr = cell_id.array(mfi);
#ifdef AMREX_USE_EB
        auto fabtyp = (flags) ? (*flags)[mfi].getType(bx) : FabType::regular;
        if (fabtyp == FabType::covered)
        {
            ncells_grid[mfi] = 0;
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D_FUSIBLE(gbx, i, j, k,
            {
                cid_arr(i,j,k) = std::numeric_limits<PetscInt>::lowest();
            });
        }
        else
#endif
        {
            Long npts = bx.numPts();
            ncells_grid[mfi] = npts;
            ncells_proc += npts;

            AMREX_HOST_DEVICE_PARALLEL_FOR_3D_FUSIBLE(gbx, i, j, k,
            {
                if (bx.contains(i,j,k)) {
                    cid_arr(i,j,k) = bx.index(IntVect{AMREX_D_DECL(i,j,k)});
                } else {
                    cid_arr(i,j,k) = std::numeric_limits<PetscInt>::lowest();
                }
            });
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

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(cell_id,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box const& bx = mfi.tilebox();
        auto os = offset[mfi];
        Array4<PetscInt> const& cid_arr = cell_id.array(mfi);
        Array4<PetscInt> const& cid_vec = cell_id_vec.array(mfi);
        AMREX_HOST_DEVICE_PARALLEL_FOR_3D_FUSIBLE(bx, i, j, k,
        {
            cid_arr(i,j,k) += os;
            cid_vec(i,j,k) = cid_arr(i,j,k);
        });
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
    const auto dx = geom.CellSizeArray();
    const int bho = (m_maxorder > 2) ? 1 : 0;
    BaseFab<PetscInt> ncols_fab;

    BaseFab<Real> mat_aos_fab, mat_vec_fab;
    BaseFab<PetscInt> cols_aos_fab, cols_vec_fab;

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
            ncols_fab.resize(bx);

            const PetscInt max_stencil_size = (fabtyp == FabType::regular) ?
                regular_stencil_size : eb_stencil_size;

            mat_aos_fab.resize(bx,max_stencil_size);
            cols_aos_fab.resize(bx,max_stencil_size);

            Array4<PetscInt> const& ncols_a = ncols_fab.array();
            Array4<PetscInt const> const& cid_a = cell_id.const_array(mfi);

            Array4<Real const> const& afab = acoefs.const_array(mfi);
            GpuArray<Array4<Real const>, AMREX_SPACEDIM> bfabs {
                AMREX_D_DECL(bcoefs[0].const_array(mfi),
                             bcoefs[1].const_array(mfi),
                             bcoefs[2].const_array(mfi))};
            Array4<Real> const& diaginvfab = diaginv.array(mfi);
            GpuArray<int,AMREX_SPACEDIM*2> bctype;
            GpuArray<Real,AMREX_SPACEDIM*2> bcl;
            for (OrientationIter oit; oit; oit++)
            {
                int cdir = oit();
                bctype[cdir] = m_bndry->bndryConds(mfi)[cdir][0];
                bcl[cdir] = m_bndry->bndryLocs(mfi)[cdir];
            }

            Real sa = scalar_a;
            Real sb = scalar_b;

            if (fabtyp == FabType::regular)
            {
                constexpr int stencil_size = 2*AMREX_SPACEDIM+1;
                BaseFab<GpuArray<Real,stencil_size> > tmpmatfab
                    (bx, 1, (GpuArray<Real,stencil_size>*)mat_aos_fab.dataPtr());

                amrex::fill(tmpmatfab,
                [=] AMREX_GPU_HOST_DEVICE (GpuArray<Real,stencil_size>& sten,
                                           int i, int j, int k)
                {
                    habec_ijmat(sten, ncols_a, diaginvfab, i, j, k, cid_a,
                                sa, afab, sb, dx, bfabs, bctype, bcl, bho,
                                Array4<int const>());
                });

                BaseFab<GpuArray<PetscInt,stencil_size> > tmpcolfab
                    (bx, 1, (GpuArray<PetscInt,stencil_size>*)cols_aos_fab.dataPtr());

                amrex::fill(tmpcolfab,
                [=] AMREX_GPU_HOST_DEVICE (GpuArray<PetscInt,stencil_size>& sten,
                                           int i, int j, int k)
                {
                    habec_cols(sten, i, j, k, cid_a);
                });
            }
#ifdef AMREX_USE_EB
            else
            {
                auto const& flag_a = flags->const_array(mfi);
                auto const& vfrac_a = vfrac->const_array(mfi);
                AMREX_D_TERM(auto const& apx = area[0]->const_array(mfi);,
                             auto const& apy = area[1]->const_array(mfi);,
                             auto const& apz = area[2]->const_array(mfi);)
                AMREX_D_TERM(auto const& fcx = fcent[0]->const_array(mfi);,
                             auto const& fcy = fcent[1]->const_array(mfi);,
                             auto const& fcz = fcent[2]->const_array(mfi);)
                auto const& barea_a = barea->const_array(mfi);
                auto const& bcent_a = bcent->const_array(mfi);
                Array4<Real const> beb = (m_eb_b_coeffs) ? m_eb_b_coeffs->const_array(mfi)
                                                         : Array4<Real const>();

                constexpr int stencil_size = AMREX_D_TERM(3,*3,*3);
                BaseFab<GpuArray<Real,stencil_size> > tmpmatfab
                    (bx, 1, (GpuArray<Real,stencil_size>*)mat_aos_fab.dataPtr());

                amrex::fill(tmpmatfab,
                [=] AMREX_GPU_HOST_DEVICE (GpuArray<Real,stencil_size>& sten,
                                           int i, int j, int k)
                {
                    habec_ijmat_eb(sten, ncols_a, diaginvfab, i, j, k, cid_a,
                                   sa, afab, sb, dx, bfabs, bctype, bcl, bho,
                                   flag_a, vfrac_a, AMREX_D_DECL(apx,apy,apz),
                                   AMREX_D_DECL(fcx,fcy,fcz),barea_a,bcent_a,beb);
                });

                BaseFab<GpuArray<PetscInt,stencil_size> > tmpcolfab
                    (bx, 1, (GpuArray<PetscInt,stencil_size>*)cols_aos_fab.dataPtr());

                amrex::fill(tmpcolfab,
                [=] AMREX_GPU_HOST_DEVICE (GpuArray<PetscInt,stencil_size>& sten,
                                           int i, int j, int k)
                {
                    habec_cols_eb(sten, i, j, k, cid_a);
                });
            }
#endif

            const PetscInt nrows = ncells_grid[mfi];
            AMREX_ASSERT(nrows == static_cast<PetscInt>(bx.numPts()));
            PetscInt* rows = cell_id_vec[mfi].dataPtr();
            PetscInt* ncols = ncols_fab.dataPtr();

            // Remove invalid elements
#if defined(AMREX_DEBUG) || defined(AMREX_TESTING)
            if (sizeof(PetscInt) < sizeof(Long)) {
                Long ntot = static_cast<Long>(nrows)*max_stencil_size;
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ntot <  static_cast<Long>(std::numeric_limits<PetscInt>::max()),
                                                 "Integer overflow: please configure Hypre with --enable-bigint");
            }
#endif
            PetscInt nelems = nrows * max_stencil_size;
            PetscInt const* cols_in = cols_aos_fab.dataPtr();
            Real const* mat_in = mat_aos_fab.dataPtr();
            PetscInt* cols;
            Real* mat;
#ifdef AMREX_USE_GPU
            if (Gpu::inLaunchRegion()) {
                cols_vec_fab.resize(bx,max_stencil_size);
                mat_vec_fab.resize(bx,max_stencil_size);
                cols = cols_vec_fab.dataPtr();
                mat = mat_vec_fab.dataPtr();
                Scan::PrefixSum<PetscInt>(nelems,
                    [=] AMREX_GPU_DEVICE (PetscInt i) -> PetscInt { return (mat_in[i] != Real(0.0)); },
                    [=] AMREX_GPU_DEVICE (PetscInt i, PetscInt const& x) {
                        if (mat_in[i] != Real(0.0)) {
                            cols[x] = cols_in[i];
                            mat[x] = mat_in[i];
                        }
                    }, Scan::Type::exclusive);
            } else
#endif
            {
                cols = const_cast<PetscInt*>(cols_in);
                mat = const_cast<Real*>(mat_in);
                PetscInt m = 0;
                for (PetscInt ielem = 0; ielem < nelems; ++ielem) {
                    if (mat_in[ielem] != Real(0.0)) {
                        if (m != ielem) {
                            cols[m] = cols_in[ielem];
                            mat[m] = mat_in[ielem];
                        }
                        ++m;
                    }
                }
            }

            Gpu::synchronize();

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

    MultiFab rhs_diag(rhs.boxArray(), rhs.DistributionMap(), 1, 0);
#ifdef AMREX_USE_OMP
#pragma omp paralle if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(rhs_diag,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& reg = mfi.validbox();
        Array4<Real> const& rhs_diag_a = rhs_diag.array(mfi);
        Array4<Real const> const& rhs_a = rhs.const_array(mfi);
        Array4<Real const> const& diaginv_a = diaginv.const_array(mfi);
#ifdef AMREX_USE_EB
        auto fabtyp = (flags) ? (*flags)[mfi].getType(reg) : FabType::regular;
        if (fabtyp == FabType::singlevalued) {
            auto const& flag = flags->const_array(mfi);
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D_FUSIBLE(reg, i, j, k,
            {
                rhs_diag_a(i,j,k) = (flag(i,j,k).isCovered()) ?
                    Real(0.0) : rhs_a(i,j,k) * diaginv_a(i,j,k);
            });
        } else if (fabtyp == FabType::regular)
#endif
        {
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D_FUSIBLE(reg, i, j, k,
            {
                rhs_diag_a(i,j,k) = rhs_a(i,j,k) * diaginv_a(i,j,k);
            });
        }
    }

    for (MFIter mfi(soln); mfi.isValid(); ++mfi)
    {
        const PetscInt nrows = ncells_grid[mfi];

        if (nrows > 0)
        {
            // soln has been set to zero.
            VecSetValues(x, nrows, cell_id_vec[mfi].dataPtr(), soln[mfi].dataPtr(), INSERT_VALUES);
            VecSetValues(b, nrows, cell_id_vec[mfi].dataPtr(), rhs_diag[mfi].dataPtr(), INSERT_VALUES);
        }
    }
}

void
PETScABecLap::getSolution (MultiFab& a_soln)
{
    bool use_tmp_mf = (a_soln.nGrowVect() != 0);
    MultiFab* l_soln = &a_soln;
    MultiFab tmp;
    if (use_tmp_mf) {
        tmp.define(a_soln.boxArray(), a_soln.DistributionMap(), 1, 0);
        l_soln = &tmp;
    }

    for (MFIter mfi(*l_soln); mfi.isValid(); ++mfi)
    {
        const PetscInt nrows = ncells_grid[mfi];
        if (nrows > 0) {
            VecGetValues(x, nrows, cell_id_vec[mfi].dataPtr(), (*l_soln)[mfi].dataPtr());
        } else {
            (*l_soln)[mfi].setVal<RunOn::Device>(0.0);
        }
    }

    if (use_tmp_mf) {
        MultiFab::Copy(a_soln, tmp, 0, 0, 1, 0);
    }
}

}
