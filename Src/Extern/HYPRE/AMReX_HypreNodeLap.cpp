#include <AMReX_HypreNodeLap.H>
#include <AMReX_VisMF.H>
#include <AMReX_MLNodeLaplacian.H>

#ifdef AMREX_USE_EB
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_MultiCutFab.H>
#include <AMReX_EBFabFactory.H>
#endif

#include <cmath>
#include <numeric>
#include <limits>
#include <type_traits>

namespace amrex {

HypreNodeLap::HypreNodeLap (const BoxArray& grids_, const DistributionMapping& dmap_,
                            const Geometry& geom_, const FabFactory<FArrayBox>& factory_,
                            const iMultiFab& owner_mask_, const iMultiFab& dirichlet_mask_,
                            MPI_Comm comm_, MLNodeLaplacian const* linop_, int verbose_)
    : grids(grids_), dmap(dmap_), geom(geom_), factory(&factory_),
      owner_mask(&owner_mask_), dirichlet_mask(&dirichlet_mask_),
      comm(comm_), linop(linop_), verbose(verbose_)
{
    Gpu::LaunchSafeGuard lsg(false); // xxxxx TODO: gpu

    static_assert(AMREX_SPACEDIM > 1, "HypreNodeLap: 1D not supported");
    static_assert(std::is_same<Real, HYPRE_Real>::value, "amrex::Real != HYPRE_Real");

    int num_procs, myid;
    MPI_Comm_size(comm, &num_procs);
    MPI_Comm_rank(comm, &myid);

    const BoxArray& nba = amrex::convert(grids,IntVect::TheNodeVector());

#if defined(AMREX_DEBUG) || defined(AMREX_TESTING)
    if (sizeof(Int) < sizeof(long)) {
        long nnodes_grids = nba.numPts();
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(nnodes_grids < static_cast<long>(std::numeric_limits<Int>::max()),
                                         "You might need to configure Hypre with --enable-bigint");
    }
#endif

    // how many non-covered nodes do we have?
    nnodes_grid.define(grids,dmap);
    node_id.define(nba,dmap,1,1);
    node_id_vec.define(grids,dmap);

    node_id.setVal(std::numeric_limits<Int>::lowest());

    Int nnodes_proc = 0;

#ifdef AMREX_USE_EB
    auto ebfactory = dynamic_cast<EBFArrayBoxFactory const*>(factory);
    if (ebfactory)
    {
        const FabArray<EBCellFlagFab>& flags = ebfactory->getMultiEBCellFlagFab();
#ifdef _OPENMP
#pragma omp parallel reduction(+:nnodes_proc)
#endif
        for (MFIter mfi(node_id); mfi.isValid(); ++mfi)
        {
            const Box& ndbx = mfi.validbox();
            const auto& nid = node_id.array(mfi);
            const auto& flag = flags.array(mfi);
            const auto& owner = owner_mask->array(mfi);
            const auto& dirichlet = dirichlet_mask->array(mfi);
            int id = 0;
            const auto lo = amrex::lbound(ndbx);
            const auto hi = amrex::ubound(ndbx);
            for         (int k = lo.z; k <= hi.z; ++k) {
                for     (int j = lo.y; j <= hi.y; ++j) {
                    for (int i = lo.x; i <= hi.x; ++i) {
                        if (!owner(i,j,k) or dirichlet(i,j,k))
                        {
                            nid(i,j,k) = std::numeric_limits<Int>::lowest();
                        }
#if (AMREX_SPACEDIM == 2)
                        else if (flag(i-1,j-1,k).isCovered() and
                                 flag(i  ,j-1,k).isCovered() and
                                 flag(i-1,j  ,k).isCovered() and
                                 flag(i  ,j  ,k).isCovered())
#endif
#if (AMREX_SPACEDIM == 3)
                        else if (flag(i-1,j-1,k-1).isCovered() and
                                 flag(i  ,j-1,k-1).isCovered() and
                                 flag(i-1,j  ,k-1).isCovered() and
                                 flag(i  ,j  ,k-1).isCovered() and
                                 flag(i-1,j-1,k  ).isCovered() and
                                 flag(i  ,j-1,k  ).isCovered() and
                                 flag(i-1,j  ,k  ).isCovered() and
                                 flag(i  ,j  ,k  ).isCovered())
#endif
                        {
                            nid(i,j,k) = std::numeric_limits<Int>::lowest();
                        }
                        else
                        {
                            nid(i,j,k) = id++;
                        }
                    }
                }
            }
            nnodes_grid[mfi] = id;
            nnodes_proc += id;
        }
    }
    else
#endif
    {
#ifdef _OPENMP
#pragma omp parallel reduction(+:nnodes_proc)
#endif
        for (MFIter mfi(node_id); mfi.isValid(); ++mfi)
        {
            const Box& ndbx = mfi.validbox();
            const auto& nid = node_id.array(mfi);
            const auto& owner = owner_mask->array(mfi);
            const auto& dirichlet = dirichlet_mask->array(mfi);
            int id = 0;
            const auto lo = amrex::lbound(ndbx);
            const auto hi = amrex::ubound(ndbx);
            for         (int k = lo.z; k <= hi.z; ++k) {
                for     (int j = lo.y; j <= hi.y; ++j) {
                    for (int i = lo.x; i <= hi.x; ++i) {
                        if (!owner(i,j,k) or dirichlet(i,j,k))
                        {
                            nid(i,j,k) = std::numeric_limits<Int>::lowest();
                        }
                        else
                        {
                            nid(i,j,k) = id++;
                        }
                    }
                }
            }
            nnodes_grid[mfi] = id;
            nnodes_proc += id;
        }
    }

    Vector<Int> nnodes_allprocs(num_procs);
    MPI_Allgather(&nnodes_proc, sizeof(Int), MPI_CHAR,
                  nnodes_allprocs.data(), sizeof(Int), MPI_CHAR,
                  comm);
    Int proc_begin = 0;
    for (int i = 0; i < myid; ++i) {
        proc_begin += nnodes_allprocs[i];
    }

#ifdef AMREX_DEBUG
    Int nnodes_total = 0;
    for (auto n : nnodes_allprocs) {
        nnodes_total += n;
    }
#endif

    LayoutData<Int> offset(grids,dmap);
    Int proc_end = proc_begin;
    for (MFIter mfi(nnodes_grid); mfi.isValid(); ++mfi)
    {
        offset[mfi] = proc_end;
        proc_end += nnodes_grid[mfi];
    }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(proc_end == proc_begin+nnodes_proc,
                                     "HypreNodeLap: how did this happen?");

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(node_id,true); mfi.isValid(); ++mfi)
    {
        node_id[mfi].plus(offset[mfi], mfi.tilebox());
    }    

    node_id.FillBoundary(geom.periodicity());

    // Create and initialize A, b & x
    Int ilower = proc_begin;
    Int iupper = proc_end-1;

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

    Vector<Int> ncols;
    Vector<Int> cols;
    Vector<Real> mat;
    constexpr Int max_stencil_size = AMREX_D_TERM(3,*3,*3);

    for (MFIter mfi(node_id); mfi.isValid(); ++mfi)
    {
        const Int nrows = nnodes_grid[mfi];
        if (nrows > 0)
        {
            ncols.clear();
            ncols.reserve(nrows);

            Vector<Int> rows = node_id_vec[mfi];
            rows.clear();
            rows.reserve(nrows);

            cols.clear();
            cols.reserve(nrows*max_stencil_size);

            mat.clear();
            mat.reserve(nrows*max_stencil_size);

            const Array4<Int const> nid = node_id.array(mfi);

            linop->fillIJMatrix(mfi, nid, ncols, rows, cols, mat);

#ifdef AMREX_DEBUG
            Int nvalues = 0;
            for (Int i = 0; i < nrows; ++i) {
                nvalues += ncols[i];
            }
            for (Int i = 0; i < nvalues; ++i) {
                AMREX_ASSERT(cols[i] >= 0 && cols[i] < nnodes_total);
            }
#endif

            HYPRE_IJMatrixSetValues(A, nrows, ncols.data(), rows.data(), cols.data(), mat.data());
        }
    }
    HYPRE_IJMatrixAssemble(A);

    // Create solver
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

    HYPRE_ParCSRMatrix par_A = NULL;
    HYPRE_IJMatrixGetObject(A, (void**)  &par_A);
    HYPRE_BoomerAMGSetup(solver, par_A, NULL, NULL);
}

HypreNodeLap::~HypreNodeLap ()
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

}
