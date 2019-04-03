#include <AMReX_HypreNodeLap.H>
#include <AMReX_VisMF.H>

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
                            MPI_Comm comm_)
    : grids(grids_), dmap(dmap_), geom(geom_), factory(&factory_), comm(comm_),
      owner_mask(&owner_mask_), dirichlet_mask(&dirichlet_mask_)
{
    Gpu::LaunchSafeGuard lsg(false); // xxxxx TODO: gpu

    static_assert(AMREX_SPACEDIM > 1, "HypreNodeLap: 1D not supported");
    static_assert(std::is_same<Real, HYPRE_Real>::value, "amrex::Real != HYPRE_Real");

    int num_procs, myid;
    MPI_Comm_size(comm, &num_procs);
    MPI_Comm_rank(comm, &myid);

    const BoxArray& nba = amrex::convert(grids,IntVect::TheNodeVector());

#if defined(AMREX_DEBUG) || defined(AMREX_TESTING)
    if (sizeof(HYPRE_Int) < sizeof(long)) {
        long nnodes_grids = nba.numPts();
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(nnodes_grids < static_cast<long>(std::numeric_limits<HYPRE_Int>::max()),
                                         "You might need to configure Hypre with --enable-bigint");
    }
#endif

    // how many non-covered nodes do we have?
    nnodes_grid.define(grids,dmap);
    node_id.define(nba,dmap,1,1);
    node_id_vec.define(grids,dmap);

    HYPRE_Int nnodes_proc = 0;

//    int i = 0;
//    amrex::Print() << "i = " << i << std::endl;

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
                            nid(i,j,k) = std::numeric_limits<HYPRE_Int>::lowest();
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
                            nid(i,j,k) = std::numeric_limits<HYPRE_Int>::lowest();
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
                            nid(i,j,k) = std::numeric_limits<HYPRE_Int>::lowest();
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

    Vector<HYPRE_Int> nnodes_allprocs(num_procs);
    MPI_Allgather(&nnodes_proc, sizeof(HYPRE_Int), MPI_CHAR,
                  nnodes_allprocs.data(), sizeof(HYPRE_Int), MPI_CHAR,
                  comm);
    HYPRE_Int proc_begin = 0;
    for (int i = 0; i < myid; ++i) {
        proc_begin += nnodes_allprocs[i];
    }

#ifdef AMREX_DEBUG
    HYPRE_Int nnodes_total = 0;
    for (auto n : nnodes_allprocs) {
        nnodes_total += n;
    }
#endif

    LayoutData<HYPRE_Int> offset(grids,dmap);
    HYPRE_Int proc_end = proc_begin;
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

    for (MFIter mfi(node_id); mfi.isValid(); ++mfi)
    {
        constexpr HYPRE_Int max_stencil_size = AMREX_D_TERM(3,*3,*3);
    }

    amrex::Abort("HypreNodeLap ctor");
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
