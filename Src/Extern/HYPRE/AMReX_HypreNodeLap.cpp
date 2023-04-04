#include <AMReX_HypreNodeLap.H>
#include <AMReX_VisMF.H>
#include <AMReX_MLNodeLaplacian.H>

#include <cmath>
#include <numeric>
#include <limits>
#include <type_traits>

namespace amrex {

HypreNodeLap::HypreNodeLap (const BoxArray& grids_, const DistributionMapping& dmap_, // NOLINT(modernize-pass-by-value)
                            const Geometry& geom_, const FabFactory<FArrayBox>& factory_,
                            const iMultiFab& owner_mask_, const iMultiFab& dirichlet_mask_,
                            MPI_Comm comm_, MLNodeLinOp const* linop_, int verbose_,
                            std::string options_namespace_)
    : grids(grids_), dmap(dmap_), geom(geom_), factory(&factory_),
      owner_mask(&owner_mask_), dirichlet_mask(&dirichlet_mask_),
      comm(comm_), linop(linop_), verbose(verbose_),
      options_namespace(std::move(options_namespace_))
{
    static_assert(AMREX_SPACEDIM > 1, "HypreNodeLap: 1D not supported");
    static_assert(std::is_same<Real, HYPRE_Real>::value, "amrex::Real != HYPRE_Real");

    int num_procs, myid;
    MPI_Comm_size(comm, &num_procs);
    MPI_Comm_rank(comm, &myid);

    const BoxArray& nba = amrex::convert(grids,IntVect::TheNodeVector());

#if defined(AMREX_DEBUG) || defined(AMREX_TESTING)
    if (sizeof(Int) < sizeof(Long)) {
        Long nnodes_grids = nba.numPts();
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(nnodes_grids < static_cast<Long>(std::numeric_limits<Int>::max()),
                                         "You might need to configure Hypre with --enable-bigint");
    }
#endif

    // how many non-covered nodes do we have?
    nnodes_grid.define(nba,dmap);
    node_id.define(nba,dmap,1,1);
    node_id_vec.define(nba,dmap);
    local_node_id.define(nba,dmap,1,0);
    id_offset.define(nba,dmap);
    tmpsoln.define(nba,dmap,1,0);

    Int nnodes_proc = fill_local_node_id();

    // At this point, local_node_id stores the ids local to each box.
    // nnodes_grid stroes the number of nodes in each box.  nnodes_proc is
    // the number of nodes on this MPI process.  If a nodal is invalid, its
    // id is invalid (i.e., a very negative number).  Note that the data
    // type of local_node_id is int, not HYPRE_Int for performance on GPU.

    Vector<Int> nnodes_allprocs(num_procs);
    MPI_Allgather(&nnodes_proc, sizeof(Int), MPI_CHAR,
                  nnodes_allprocs.data(), sizeof(Int), MPI_CHAR,
                  comm);
    Int proc_begin = 0;
    for (int i = 0; i < myid; ++i) {
        proc_begin += nnodes_allprocs[i];
    }

    Int proc_end = proc_begin;
    for (MFIter mfi(nnodes_grid); mfi.isValid(); ++mfi)
    {
        id_offset[mfi] = proc_end;
        proc_end += nnodes_grid[mfi];
    }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(proc_end == proc_begin+nnodes_proc,
                                     "HypreNodeLap: how did this happen?");

    fill_global_node_id(); // node_id = local_node_id + id_offset, and fill node_id_vec

    amrex::OverrideSync(node_id, *owner_mask, geom.periodicity());
    node_id.FillBoundary(geom.periodicity());

    // Create and initialize A, b & x
    Int ilower = proc_begin;
    Int iupper = proc_end-1;

    hypre_ij = std::make_unique<HypreIJIface>(comm, ilower, iupper, verbose);
    hypre_ij->parse_inputs(options_namespace);

    // Obtain non-owning references to the matrix, rhs, and solution data
    A = hypre_ij->A();
    b = hypre_ij->b();
    x = hypre_ij->x();

    Gpu::DeviceVector<Int> ncols_vec;
    Gpu::DeviceVector<Int> cols_vec;
    Gpu::DeviceVector<Real> mat_vec;
    constexpr Int max_stencil_size = AMREX_D_TERM(3,*3,*3);

    for (MFIter mfi(node_id, MFItInfo{}.UseDefaultStream()); mfi.isValid(); ++mfi)
    {
        const Int nrows = nnodes_grid[mfi];
        if (nrows > 0)
        {
            ncols_vec.clear();
            ncols_vec.resize(nrows);
            auto* ncols = ncols_vec.data();

            const auto& rows_vec = node_id_vec[mfi];
            const auto* rows = rows_vec.data();

            cols_vec.clear();
            cols_vec.resize(nrows*max_stencil_size);
            auto* cols = cols_vec.data();

            mat_vec.clear();
            mat_vec.resize(nrows*max_stencil_size);
            auto* mat = mat_vec.data();

            const auto& gid = node_id.const_array(mfi);
            const auto& lid = local_node_id.const_array(mfi);

            linop->fillIJMatrix(mfi, gid, lid, ncols, cols, mat);

            if (hypre_ij->adjustSingularMatrix() && linop->isBottomSingular()
                && id_offset[mfi] == 0 && nnodes_grid[mfi] > 0)
            {
                adjust_singular_matrix(ncols, cols, rows, mat);
            }

            Gpu::synchronize();
            HYPRE_IJMatrixSetValues(A, nrows, ncols, rows, cols, mat);
            Gpu::synchronize();
        }
    }
    HYPRE_IJMatrixAssemble(A);
}

HypreNodeLap::~HypreNodeLap () = default;

void
HypreNodeLap::solve (MultiFab& soln, const MultiFab& rhs,
                     Real rel_tol, Real abs_tol, int max_iter)
{
    BL_PROFILE("HypreNodeLap::solve()");

    HYPRE_IJVectorInitialize(b);
    HYPRE_IJVectorInitialize(x);
    //
    loadVectors(soln, rhs);
    //
    HYPRE_IJVectorAssemble(x);
    HYPRE_IJVectorAssemble(b);

    hypre_ij->solve(rel_tol, abs_tol, max_iter);

    getSolution(soln);
}

HypreNodeLap::Int
HypreNodeLap::fill_local_node_id ()
{
#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion()) {
        return fill_local_node_id_gpu();
    } else
#endif
    {
        return fill_local_node_id_cpu();
    }
}

#ifdef AMREX_USE_GPU

HypreNodeLap::Int
HypreNodeLap::fill_local_node_id_gpu ()
{
    Int nnodes_proc = 0;

    for (MFIter mfi(local_node_id); mfi.isValid(); ++mfi)
    {
        const Box& ndbx = mfi.validbox();
        const auto& nid = local_node_id.array(mfi);
        const auto& owner = owner_mask->const_array(mfi);
        const auto& dirichlet = dirichlet_mask->const_array(mfi);
        AMREX_ASSERT(ndbx.numPts() < static_cast<Long>(std::numeric_limits<int>::max()));
        const int npts = ndbx.numPts();
        int nnodes_box = amrex::Scan::PrefixSum<int>(npts,
            [=] AMREX_GPU_DEVICE (int offset) noexcept
            {
                int valid_node = 1;
                const Dim3 cell = ndbx.atOffset(offset).dim3();
                const int i = cell.x;
                const int j = cell.y;
                const int k = cell.z;
                if (!owner(i,j,k) || dirichlet(i,j,k)) {
                    valid_node = 0;
                }
                nid(i,j,k) = valid_node;
                return valid_node;
            },
            [=] AMREX_GPU_DEVICE (int offset, int ps) noexcept
            {
                const Dim3 cell = ndbx.atOffset(offset).dim3();
                nid(cell.x,cell.y,cell.z) = nid(cell.x,cell.y,cell.z)
                    ? ps : std::numeric_limits<int>::lowest();
            },
            amrex::Scan::Type::exclusive);

        nnodes_grid[mfi] = nnodes_box;
        nnodes_proc += nnodes_box;
    }

    return nnodes_proc;
}

#endif

HypreNodeLap::Int
HypreNodeLap::fill_local_node_id_cpu ()
{
    Int nnodes_proc = 0;

#ifdef AMREX_USE_OMP
#pragma omp parallel reduction(+:nnodes_proc)
#endif
    for (MFIter mfi(local_node_id); mfi.isValid(); ++mfi)
    {
        const Box& ndbx = mfi.validbox();
        const auto& nid = local_node_id.array(mfi);
        const auto& owner = owner_mask->const_array(mfi);
        const auto& dirichlet = dirichlet_mask->const_array(mfi);
        int id = 0;
        const auto lo = amrex::lbound(ndbx);
        const auto hi = amrex::ubound(ndbx);
        for         (int k = lo.z; k <= hi.z; ++k) {
            for     (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    if (!owner(i,j,k) || dirichlet(i,j,k))
                    {
                        nid(i,j,k) = std::numeric_limits<int>::lowest();
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

    return nnodes_proc;
}

void
HypreNodeLap::fill_global_node_id ()
{
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(node_id); mfi.isValid(); ++mfi)
    {
        Int os = id_offset[mfi];
        const Box& bx = mfi.growntilebox();
        const auto& nid = node_id.array(mfi);
        const auto& lnid = local_node_id.const_array(mfi);
        auto& rows_vec = node_id_vec[mfi];
        rows_vec.resize(nnodes_grid[mfi]);
        auto* rows = rows_vec.data();
        AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
        {
            if (lnid.contains(i,j,k) && lnid(i,j,k) >= 0) {
                const auto gid = lnid(i,j,k) + os;
                rows[lnid(i,j,k)] = gid;
                nid(i,j,k) = static_cast<AtomicInt>(gid);
            } else {
                nid(i,j,k) = std::numeric_limits<AtomicInt>::max();
            }
        });
    }
}

void
HypreNodeLap::adjust_singular_matrix (Int const* ncols, Int const* cols, Int const* rows, Real* mat)
{
    AMREX_HOST_DEVICE_FOR_1D(1, m,
    {
        amrex::ignore_unused(m);
        const int num_cols = ncols[0];
        for (int ic = 0; ic < num_cols; ++ic) {
            mat[ic] = (cols[ic] == rows[0]) ? mat[ic] : Real(0.0);
        }
    });
}

void
HypreNodeLap::loadVectors (MultiFab& soln, const MultiFab& rhs)
{
    BL_PROFILE("HypreNodeLap::loadVectors()");

    soln.setVal(0.0);

    Gpu::DeviceVector<Real> bvec;
    for (MFIter mfi(soln, MFItInfo{}.UseDefaultStream()); mfi.isValid(); ++mfi)
    {
        const Int nrows = nnodes_grid[mfi];
        if (nrows >= 0)
        {
            const auto& rows_vec = node_id_vec[mfi];
            HYPRE_IJVectorSetValues(x, nrows, rows_vec.data(), soln[mfi].dataPtr());

            bvec.clear();
            bvec.resize(nrows);
            auto* bp = bvec.data();

            const auto& bfab = rhs.array(mfi);
            const auto& lid = local_node_id.array(mfi);
            linop->fillRHS(mfi, lid, bp, bfab);

            if (hypre_ij->adjustSingularMatrix() && linop->isBottomSingular()
                && id_offset[mfi] == 0 && nnodes_grid[mfi] > 0)
            {
                AMREX_HOST_DEVICE_FOR_1D(1, m,
                {
                    amrex::ignore_unused(m);
                    bp[0] = 0.0;
                });
            }

            Gpu::synchronize();
            HYPRE_IJVectorSetValues(b, nrows, rows_vec.data(), bvec.data());
            Gpu::synchronize();
        }
    }
}

void
HypreNodeLap::getSolution (MultiFab& soln)
{
    tmpsoln.setVal(0.0);

    Gpu::DeviceVector<Real> xvec;
    for (MFIter mfi(tmpsoln, MFItInfo{}.UseDefaultStream()); mfi.isValid(); ++mfi)
    {
        const Int nrows = nnodes_grid[mfi];
        if (nrows >= 0)
        {
            const auto& rows_vec = node_id_vec[mfi];
            xvec.clear();
            xvec.resize(nrows);
            Real* xp = xvec.data();
            HYPRE_IJVectorGetValues(x, nrows, rows_vec.data(), xp);
            Gpu::synchronize();

            const Box& bx = mfi.validbox();
            const auto& xfab = tmpsoln.array(mfi);
            const auto& lid = local_node_id.array(mfi);
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
            {
                if (lid(i,j,k,0) >= 0) {
                    xfab(i,j,k) = xp[lid(i,j,k)];
                }
            });

            Gpu::synchronize();
        }
    }

    soln.ParallelAdd(tmpsoln, 0, 0, 1, geom.periodicity());
}

}
