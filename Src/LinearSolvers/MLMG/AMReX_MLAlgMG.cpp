#include <AMReX_MLAlgMG.H>
#include <AMReX_AlgMG.H>
#include <AMReX_AlgVector.H>
#include <AMReX_SpMatrix.H>
#include <AMReX_MLNodeLinOp.H>
#include <AMReX_Habec_K.H>
#if defined(AMREX_USE_EB) && (AMREX_SPACEDIM > 1)
#include <AMReX_EBFabFactory.H>
#include <AMReX_EBMultiFabUtil.H>
#endif
#include <AMReX_LayoutData.H>
#include <AMReX_ParallelReduce.H>
#include <AMReX_Reduce.H>
#include <AMReX_Scan.H>
#include <AMReX_MultiFabUtil.H>

#include <limits>
#include <numeric>

namespace amrex {

namespace {
    // A periodic direction two cells wide at this level gives a row two
    // entries with the same column. Sum them and invalidate the second.
    bool hasDuplicateColumns (Geometry const& geom)
    {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            if (geom.isPeriodic(idim) && geom.Domain().length(idim) <= 2) { return true; }
        }
        return false;
    }

    void mergeDuplicateColumns (Gpu::DeviceVector<Real>& mat, Gpu::DeviceVector<Long>& cols,
                                Gpu::DeviceVector<Long> const& row_offset, Long nrows)
    {
        Real* AMREX_RESTRICT pm = mat.data();
        Long* AMREX_RESTRICT pc = cols.data();
        Long const* AMREX_RESTRICT ro = row_offset.data();
        ParallelFor(nrows, [=] AMREX_GPU_DEVICE (Long r) noexcept
        {
            for (Long a = ro[r]; a < ro[r+1]; ++a) {
                if (pc[a] < 0) { continue; }
                for (Long b = a+1; b < ro[r+1]; ++b) {
                    if (pc[b] == pc[a]) {
                        pm[a] += pm[b];
                        pm[b] = Real(0.0);
                        pc[b] = -1;
                    }
                }
            }
        });
        Gpu::streamSynchronize();
    }
}

namespace {

// Rank-contiguous row numbering from per-rank counts.
AlgPartition make_partition (Long nrows_proc)
{
    int const nprocs = ParallelContext::NProcsSub();
    Vector<Long> counts(nprocs, nrows_proc);
#ifdef AMREX_USE_MPI
    ParallelAllGather::AllGather(nrows_proc, counts.data(), ParallelContext::CommunicatorSub());
#endif
    Vector<Long> rows(nprocs+1, 0); // exclusive prefix sum of counts
    std::partial_sum(counts.begin(), counts.end(), rows.begin()+1);
    return AlgPartition(std::move(rows));
}

}

struct MLAlgMG::Impl
{
    Impl (int mglev, BoxArray const& grids, DistributionMapping const& dmap,
          Geometry const& geom, iMultiFab const& owner_mask,
          iMultiFab const& dirichlet_mask, MLNodeLinOp const& linop);

    Impl (int mglev, BoxArray const& grids, DistributionMapping const& dmap,
          Geometry const& geom, FabFactory<FArrayBox> const& factory,
          iMultiFab const* overset_mask, Real ascalar, Real bscalar,
          MultiFab const* acoef, Array<MultiFab const*,AMREX_SPACEDIM> const& bcoef,
          MultiFab const* eb_bcoef,
          LayoutData<GpuArray<int,2*AMREX_SPACEDIM>> const& bctype,
          LayoutData<GpuArray<Real,2*AMREX_SPACEDIM>> const& bcl, int maxorder);

    // Constructors cannot hold device lambdas (nvcc), so they call these.
    void defineNodal (BoxArray const& grids, DistributionMapping const& dmap,
                      iMultiFab const& owner_mask, iMultiFab const& dirichlet_mask);
    void defineCell (BoxArray const& grids, DistributionMapping const& dmap,
                     FabFactory<FArrayBox> const& factory, iMultiFab const* overset_mask,
                     Real ascalar, Real bscalar, MultiFab const* acoef,
                     Array<MultiFab const*,AMREX_SPACEDIM> const& bcoef, MultiFab const* eb_bcoef,
                     LayoutData<GpuArray<int,2*AMREX_SPACEDIM>> const& bctype,
                     LayoutData<GpuArray<Real,2*AMREX_SPACEDIM>> const& bcl, int maxorder);
    void assembleNodal (MLNodeLinOp const& linop);
    void assembleCell (FabFactory<FArrayBox> const& factory, iMultiFab const* overset_mask,
                       Real ascalar, Real bscalar, MultiFab const& acoef,
                       Array<MultiFab const*,AMREX_SPACEDIM> const& bcoef,
                       MultiFab const* eb_bcoef,
                       LayoutData<GpuArray<int,2*AMREX_SPACEDIM>> const& bctype,
                       LayoutData<GpuArray<Real,2*AMREX_SPACEDIM>> const& bcl, int maxorder);
    void loadRHS (MultiFab const& rhs);
    void getSolution (MultiFab& soln);

    int m_mglev = 0;
    bool m_nodal = true;
    MLNodeLinOp const* m_nodelinop = nullptr;
    Geometry m_geom;
    iMultiFab m_lid;                     // nodal: local row id, negative for no row
    FabArray<BaseFab<Long>> m_gid;       // global row id, 1 ghost; nodal: max() for
                                         // no row, cell: lowest()
    LayoutData<Long> m_nrows_grid;
    LayoutData<Long> m_row_begin;        // first local row of each box
    Long m_nrows_proc = 0;
    MultiFab m_tmp;                      // scratch for the nodal scatter

    // cell-centered
    iMultiFab const* m_overset_mask = nullptr;
    iMultiFab m_osm_grown;               // overset mask with filled ghost cells
#if defined(AMREX_USE_EB) && (AMREX_SPACEDIM > 1)
    FabArray<EBCellFlagFab> const* m_flags = nullptr;
    MultiFab const* m_vfrac = nullptr;
#endif

    AlgPartition m_part;
    SpMatrix<Real> m_A;
    AlgVector<Real> m_x;
    AlgVector<Real> m_b;
    AlgMG<Real> m_solver;
};

MLAlgMG::MLAlgMG (int mglev, BoxArray const& grids, DistributionMapping const& dmap,
                  Geometry const& geom, iMultiFab const& owner_mask,
                  iMultiFab const& dirichlet_mask, MLNodeLinOp const& linop)
    : m_impl(std::make_unique<Impl>(mglev, grids, dmap, geom, owner_mask,
                                    dirichlet_mask, linop))
{}

MLAlgMG::~MLAlgMG () = default;

AlgMG<Real>& MLAlgMG::solver () noexcept { return m_impl->m_solver; }

void
MLAlgMG::applyVcycle (MultiFab& soln, MultiFab const& rhs)
{
    BL_PROFILE("MLAlgMG::applyVcycle()");

    m_impl->loadRHS(rhs);
    m_impl->m_solver.precond(m_impl->m_x, m_impl->m_b);
    m_impl->getSolution(soln);
}

void
MLAlgMG::solve (MultiFab& soln, MultiFab const& rhs, Real reltol, Real abstol, int maxiter)
{
    BL_PROFILE("MLAlgMG::solve()");

    m_impl->m_solver.setRelTol(reltol);
    m_impl->m_solver.setAbsTol(abstol);
    m_impl->m_solver.setMaxIter(maxiter);

    m_impl->loadRHS(rhs);
    m_impl->m_x.setVal(Real(0.0));
    m_impl->m_solver.solve(m_impl->m_x, m_impl->m_b);
    m_impl->getSolution(soln);
}

MLAlgMG::Impl::Impl (int mglev, BoxArray const& grids, DistributionMapping const& dmap,
                     Geometry const& geom, iMultiFab const& owner_mask,
                     iMultiFab const& dirichlet_mask, MLNodeLinOp const& linop)
    : m_mglev(mglev), m_nodelinop(&linop), m_geom(geom)
{
    defineNodal(grids, dmap, owner_mask, dirichlet_mask);
}

void
MLAlgMG::Impl::defineNodal (BoxArray const& grids, DistributionMapping const& dmap,
                            iMultiFab const& owner_mask, iMultiFab const& dirichlet_mask)
{
    BL_PROFILE("MLAlgMG::defineNodal()");

    const BoxArray& nba = amrex::convert(grids, IntVect::TheNodeVector());
    m_lid.define(nba, dmap, 1, 0);
    m_gid.define(nba, dmap, 1, 1);
    m_nrows_grid.define(nba, dmap);
    m_row_begin.define(nba, dmap);
    m_tmp.define(nba, dmap, 1, 0);

    // Local ids: owned, non-Dirichlet nodes in lexicographic order per box.
#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion()) {
        for (MFIter mfi(m_lid); mfi.isValid(); ++mfi) {
            const Box& ndbx = mfi.validbox();
            auto const& nid = m_lid.array(mfi);
            auto const& owner = owner_mask.const_array(mfi);
            auto const& dirichlet = dirichlet_mask.const_array(mfi);
            AMREX_ASSERT(ndbx.numPts() < static_cast<Long>(std::numeric_limits<int>::max()));
            const int npts = static_cast<int>(ndbx.numPts());
            int nnodes_box = Scan::PrefixSum<int>(npts,
                [=] AMREX_GPU_DEVICE (int offset) noexcept -> int
                {
                    const Dim3 cell = ndbx.atOffset(offset).dim3();
                    int valid = (owner(cell.x,cell.y,cell.z) && !dirichlet(cell.x,cell.y,cell.z)) ? 1 : 0;
                    nid(cell.x,cell.y,cell.z) = valid;
                    return valid;
                },
                [=] AMREX_GPU_DEVICE (int offset, int ps) noexcept
                {
                    const Dim3 cell = ndbx.atOffset(offset).dim3();
                    nid(cell.x,cell.y,cell.z) = nid(cell.x,cell.y,cell.z)
                        ? ps : std::numeric_limits<int>::lowest();
                },
                Scan::Type::exclusive);
            m_nrows_grid[mfi] = nnodes_box;
            m_nrows_proc += nnodes_box;
        }
    } else
#endif
    {
        for (MFIter mfi(m_lid); mfi.isValid(); ++mfi) {
            const Box& ndbx = mfi.validbox();
            auto const& nid = m_lid.array(mfi);
            auto const& owner = owner_mask.const_array(mfi);
            auto const& dirichlet = dirichlet_mask.const_array(mfi);
            int id = 0;
            amrex::LoopOnCpu(ndbx, [&] (int i, int j, int k) noexcept
            {
                if (!owner(i,j,k) || dirichlet(i,j,k)) {
                    nid(i,j,k) = std::numeric_limits<int>::lowest();
                } else {
                    nid(i,j,k) = id++;
                }
            });
            m_nrows_grid[mfi] = id;
            m_nrows_proc += id;
        }
    }

    m_part = make_partition(m_nrows_proc);
    Long const proc_begin = m_part.globalRowBegin();

    // Global ids, then let the owners overwrite the shared nodes.
    Long os = proc_begin;
    for (MFIter mfi(m_gid); mfi.isValid(); ++mfi) {
        m_row_begin[mfi] = os - proc_begin;
        const Box& bx = mfi.growntilebox();
        auto const& gid = m_gid.array(mfi);
        auto const& lid = m_lid.const_array(mfi);
        const Long os_box = os;
        AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
        {
            if (lid.contains(i,j,k) && lid(i,j,k) >= 0) {
                gid(i,j,k) = lid(i,j,k) + os_box;
            } else {
                gid(i,j,k) = std::numeric_limits<Long>::max();
            }
        });
        os += m_nrows_grid[mfi];
    }
    AMREX_ALWAYS_ASSERT(os == proc_begin + m_nrows_proc);

    amrex::OverrideSync(m_gid, owner_mask, m_geom.periodicity());
    m_gid.FillBoundary(m_geom.periodicity());

    assembleNodal(*m_nodelinop);

    m_x.define(m_part);
    m_b.define(m_part);
    m_solver.define(m_A);
}

void
MLAlgMG::Impl::assembleNodal (MLNodeLinOp const& linop)
{
    BL_PROFILE("MLAlgMG::assembleNodal()");

    constexpr int max_stencil = AMREX_D_TERM(3,*3,*3);

    Gpu::DeviceVector<Long> ncols(m_nrows_proc);
    Gpu::DeviceVector<Long> cols(m_nrows_proc*max_stencil);
    Gpu::DeviceVector<Real> mat(m_nrows_proc*max_stencil);
    Gpu::DeviceVector<Long> cols_box;
    Gpu::DeviceVector<Real> mat_box;

    Long nnz = 0;
    for (MFIter mfi(m_lid, MFItInfo{}.UseDefaultStream()); mfi.isValid(); ++mfi) {
        const Long nrows = m_nrows_grid[mfi];
        if (nrows == 0) { continue; }

        cols_box.clear();
        cols_box.resize(nrows*max_stencil);
        mat_box.clear();
        mat_box.resize(nrows*max_stencil);

        Long* ncols_p = ncols.data() + m_row_begin[mfi];
        linop.fillAlgMatrix(m_mglev, mfi, m_gid.const_array(mfi), m_lid.const_array(mfi),
                            ncols_p, cols_box.data(), mat_box.data());

        Long nnz_box = Reduce::Sum<Long>(nrows,
            [=] AMREX_GPU_DEVICE (Long i) -> Long { return ncols_p[i]; });
        Gpu::copyAsync(Gpu::deviceToDevice, cols_box.begin(), cols_box.begin()+nnz_box,
                       cols.begin()+nnz);
        Gpu::copyAsync(Gpu::deviceToDevice, mat_box.begin(), mat_box.begin()+nnz_box,
                       mat.begin()+nnz);
        Gpu::streamSynchronize();
        nnz += nnz_box;
    }

    Gpu::DeviceVector<Long> row_offset(m_nrows_proc+1);
    Long const total = Scan::ExclusiveSum(m_nrows_proc, ncols.data(), row_offset.data());
    Gpu::streamSynchronize();
    AMREX_ALWAYS_ASSERT(total == nnz);
    Long* last = row_offset.data() + m_nrows_proc;
    AMREX_HOST_DEVICE_FOR_1D(1, i, { amrex::ignore_unused(i); *last = total; });
    Gpu::streamSynchronize();

    bool const merged = hasDuplicateColumns(m_geom);
    if (merged) {
        mergeDuplicateColumns(mat, cols, row_offset, m_nrows_proc);
    }

    m_A.define(m_part, mat.data(), cols.data(), nnz, row_offset.data(),
               CsrSorted{false}, CsrValid{!merged});
}

MLAlgMG::MLAlgMG (int mglev, BoxArray const& grids, DistributionMapping const& dmap,
                  Geometry const& geom, FabFactory<FArrayBox> const& factory,
                  iMultiFab const* overset_mask, Real ascalar, Real bscalar,
                  MultiFab const* acoef, Array<MultiFab const*,AMREX_SPACEDIM> const& bcoef,
                  MultiFab const* eb_bcoef,
                  LayoutData<GpuArray<int,2*AMREX_SPACEDIM>> const& bctype,
                  LayoutData<GpuArray<Real,2*AMREX_SPACEDIM>> const& bcl, int maxorder)
    : m_impl(std::make_unique<Impl>(mglev, grids, dmap, geom, factory, overset_mask,
                                    ascalar, bscalar, acoef, bcoef, eb_bcoef,
                                    bctype, bcl, maxorder))
{}

MLAlgMG::Impl::Impl (int mglev, BoxArray const& grids, DistributionMapping const& dmap,
                     Geometry const& geom, FabFactory<FArrayBox> const& factory,
                     iMultiFab const* overset_mask, Real ascalar, Real bscalar,
                     MultiFab const* acoef, Array<MultiFab const*,AMREX_SPACEDIM> const& bcoef,
                     MultiFab const* eb_bcoef,
                     LayoutData<GpuArray<int,2*AMREX_SPACEDIM>> const& bctype,
                     LayoutData<GpuArray<Real,2*AMREX_SPACEDIM>> const& bcl, int maxorder)
    : m_mglev(mglev), m_nodal(false), m_geom(geom), m_overset_mask(overset_mask)
{
    defineCell(grids, dmap, factory, overset_mask, ascalar, bscalar, acoef, bcoef, eb_bcoef,
               bctype, bcl, maxorder);
}

void
MLAlgMG::Impl::defineCell (BoxArray const& grids, DistributionMapping const& dmap,
                           FabFactory<FArrayBox> const& factory, iMultiFab const* overset_mask,
                           Real ascalar, Real bscalar, MultiFab const* acoef,
                           Array<MultiFab const*,AMREX_SPACEDIM> const& bcoef, MultiFab const* eb_bcoef,
                           LayoutData<GpuArray<int,2*AMREX_SPACEDIM>> const& bctype,
                           LayoutData<GpuArray<Real,2*AMREX_SPACEDIM>> const& bcl, int maxorder)
{
    BL_PROFILE("MLAlgMG::defineCell()");

    m_gid.define(grids, dmap, 1, 1);
    m_nrows_grid.define(grids, dmap);
    m_row_begin.define(grids, dmap);

#if defined(AMREX_USE_EB) && (AMREX_SPACEDIM > 1)
    auto const* ebfactory = dynamic_cast<EBFArrayBoxFactory const*>(&factory);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(overset_mask == nullptr || ebfactory == nullptr,
                                     "MLAlgMG: cannot have both EB and overset");
    m_flags = ebfactory ? &(ebfactory->getMultiEBCellFlagFab()) : nullptr;
    m_vfrac = ebfactory ? &(ebfactory->getVolFrac()) : nullptr;
#endif

    // Rows: every cell of every box, except boxes that are fully covered.
    // Ghost cells outside the domain and cells of covered boxes get
    // lowest(), which the kernels read as "no row".
    for (MFIter mfi(m_gid); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        const Box& gbx = amrex::grow(bx,1);
        auto const& gid = m_gid.array(mfi);
#if defined(AMREX_USE_EB) && (AMREX_SPACEDIM > 1)
        auto fabtyp = m_flags ? (*m_flags)[mfi].getType(bx) : FabType::regular;
#else
        auto fabtyp = FabType::regular;
#endif
        Long const nrows = (fabtyp == FabType::covered) ? 0 : bx.numPts();
        m_nrows_grid[mfi] = nrows;
        m_nrows_proc += nrows;
        AMREX_HOST_DEVICE_PARALLEL_FOR_3D(gbx, i, j, k,
        {
            if (nrows > 0 && bx.contains(i,j,k)) {
                gid(i,j,k) = bx.index(IntVect{AMREX_D_DECL(i,j,k)});
            } else {
                gid(i,j,k) = std::numeric_limits<Long>::lowest();
            }
        });
    }

    m_part = make_partition(m_nrows_proc);
    Long const proc_begin = m_part.globalRowBegin();

    Long os = proc_begin;
    for (MFIter mfi(m_gid); mfi.isValid(); ++mfi) {
        m_row_begin[mfi] = os - proc_begin;
        if (m_nrows_grid[mfi] > 0) {
            const Box& bx = mfi.validbox();
            auto const& gid = m_gid.array(mfi);
            const Long os_box = os;
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
            {
                gid(i,j,k) += os_box;
            });
        }
        os += m_nrows_grid[mfi];
    }
    AMREX_ALWAYS_ASSERT(os == proc_begin + m_nrows_proc);

    m_gid.FillBoundary(m_geom.periodicity());

    // The known cells hold a zero correction, so the couplings into them
    // are dropped below to keep the matrix symmetric. That needs the mask
    // of the neighbors: ghost cells outside the domain count as unknown.
    if (overset_mask) {
        m_osm_grown.define(grids, dmap, 1, 1);
        m_osm_grown.setVal(1);
        iMultiFab::Copy(m_osm_grown, *overset_mask, 0, 0, 1, 0);
        m_osm_grown.FillBoundary(m_geom.periodicity());
    }

    // Default coefficients: a = 0, b = 1.
    MultiFab alpha;
    if (acoef == nullptr) {
        alpha.define(grids, dmap, 1, 0, MFInfo().SetArena(The_Async_Arena()), factory);
        alpha.setVal(Real(0.0));
        acoef = &alpha;
    }
    // The operator supplies the b coefficients (with metric terms if any).
    AMREX_ALWAYS_ASSERT(bcoef[0] != nullptr);

    assembleCell(factory, overset_mask, ascalar, bscalar, *acoef, bcoef, eb_bcoef,
                 bctype, bcl, maxorder);

    m_x.define(m_part);
    m_b.define(m_part);
    m_solver.define(m_A);
}

void
MLAlgMG::Impl::assembleCell (FabFactory<FArrayBox> const& factory, iMultiFab const* overset_mask,
                             Real ascalar, Real bscalar, MultiFab const& acoef,
                             Array<MultiFab const*,AMREX_SPACEDIM> const& bcoef,
                             MultiFab const* eb_bcoef,
                             LayoutData<GpuArray<int,2*AMREX_SPACEDIM>> const& bctype,
                             LayoutData<GpuArray<Real,2*AMREX_SPACEDIM>> const& bcl, int maxorder)
{
    BL_PROFILE("MLAlgMG::assembleCell()");

    amrex::ignore_unused(factory, eb_bcoef);

    constexpr int NR = 2*AMREX_SPACEDIM+1;
    int const reg_stencil = (maxorder > 3) ? NR+AMREX_SPACEDIM : NR;

#if defined(AMREX_USE_EB) && (AMREX_SPACEDIM > 1)
    constexpr int eb_stencil = AMREX_D_TERM(3,*3,*3);
    auto const* ebfactory = dynamic_cast<EBFArrayBoxFactory const*>(&factory);
    const MultiFab* vfrac = ebfactory ? &(ebfactory->getVolFrac()) : nullptr;
    auto area = ebfactory ? ebfactory->getAreaFrac()
        : Array<const MultiCutFab*,AMREX_SPACEDIM>{AMREX_D_DECL(nullptr,nullptr,nullptr)};
    auto fcent = ebfactory ? ebfactory->getFaceCent()
        : Array<const MultiCutFab*,AMREX_SPACEDIM>{AMREX_D_DECL(nullptr,nullptr,nullptr)};
    auto const* barea = ebfactory ? &(ebfactory->getBndryArea()) : nullptr;
    auto const* bcent = ebfactory ? &(ebfactory->getBndryCent()) : nullptr;
#endif

    // Padded stencils: each box contributes nrows*stencil entries; entries
    // with a negative column or a zero value are dropped by SpMatrix::define.
    Long nentries = 0;
    LayoutData<Long> entry_begin(m_gid.boxArray(), m_gid.DistributionMap());
    LayoutData<int> stencil_size(m_gid.boxArray(), m_gid.DistributionMap());
    for (MFIter mfi(m_gid); mfi.isValid(); ++mfi) {
        entry_begin[mfi] = nentries;
        int ss = reg_stencil;
#if defined(AMREX_USE_EB) && (AMREX_SPACEDIM > 1)
        if (m_flags && (*m_flags)[mfi].getType(mfi.validbox()) == FabType::singlevalued) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(maxorder <= 3, "MLAlgMG: EB supports maxorder <= 3");
            ss = eb_stencil;
        }
#endif
        stencil_size[mfi] = ss;
        nentries += m_nrows_grid[mfi] * ss;
    }

    Gpu::DeviceVector<Real> mat(nentries);
    Gpu::DeviceVector<Long> cols(nentries);
    Gpu::DeviceVector<Long> row_offset(m_nrows_proc+1);

    const auto dx = m_geom.CellSizeArray();

    for (MFIter mfi(m_gid); mfi.isValid(); ++mfi) {
        const Long nrows = m_nrows_grid[mfi];
        if (nrows == 0) { continue; }
        const Box& bx = mfi.validbox();
        const int ss = stencil_size[mfi];

        Long* ro = row_offset.data() + m_row_begin[mfi];
        const Long eb0 = entry_begin[mfi];
        ParallelFor(nrows, [=] AMREX_GPU_DEVICE (Long r) noexcept { ro[r] = eb0 + r*ss; });

        // The kernels count the columns per row; the counts are not needed
        // with padded stencils.
        BaseFab<Long> ncols_fab(bx, 1, The_Async_Arena());
        Array4<Long> const& ncols_a = ncols_fab.array();
        Array4<Long const> const& cid_a = m_gid.const_array(mfi);
        Array4<Real const> const& afab = acoef.const_array(mfi);
        GpuArray<Array4<Real const>, AMREX_SPACEDIM> bfabs {
            AMREX_D_DECL(bcoef[0]->const_array(mfi),
                         bcoef[1]->const_array(mfi),
                         bcoef[2]->const_array(mfi))};
        GpuArray<int,AMREX_SPACEDIM*2> const bct = bctype[mfi];
        GpuArray<Real,AMREX_SPACEDIM*2> const bcloc = bcl[mfi];
        Real const sa = ascalar;
        Real const sb = bscalar;

        Real* matp = mat.data() + entry_begin[mfi];
        Long* colp = cols.data() + entry_begin[mfi];

        if (ss == reg_stencil)
        {
            auto osmsk = overset_mask ? m_osm_grown.const_array(mfi) : Array4<int const>();
            if (maxorder > 3) {
                habec_ij_fill<NR+AMREX_SPACEDIM>(bx, matp, colp, ncols_a, cid_a, sa, afab, sb, dx, bfabs,
                                      bct, bcloc, maxorder, osmsk, true);
            } else {
                habec_ij_fill<NR>(bx, matp, colp, ncols_a, cid_a, sa, afab, sb, dx, bfabs,
                                  bct, bcloc, maxorder, osmsk, true);
            }
        }
#if defined(AMREX_USE_EB) && (AMREX_SPACEDIM > 1)
        else
        {
            auto const& flag_a = m_flags->const_array(mfi);
            auto const& vfrac_a = vfrac->const_array(mfi);
            AMREX_D_TERM(auto const& apx = area[0]->const_array(mfi);,
                         auto const& apy = area[1]->const_array(mfi);,
                         auto const& apz = area[2]->const_array(mfi);)
            AMREX_D_TERM(auto const& fcx = fcent[0]->const_array(mfi);,
                         auto const& fcy = fcent[1]->const_array(mfi);,
                         auto const& fcz = fcent[2]->const_array(mfi);)
            auto const& barea_a = barea->const_array(mfi);
            auto const& bcent_a = bcent->const_array(mfi);
            Array4<Real const> beb = eb_bcoef ? eb_bcoef->const_array(mfi) : Array4<Real const>();
            int const bho = (maxorder > 2) ? 1 : 0;

            BaseFab<GpuArray<Real,eb_stencil>> tmpmatfab
                (bx, 1, (GpuArray<Real,eb_stencil>*)matp);
            amrex::fill(tmpmatfab,
            [=] AMREX_GPU_HOST_DEVICE (GpuArray<Real,eb_stencil>& sten, int i, int j, int k)
            {
                habec_ijmat_eb(sten, ncols_a, i, j, k, bx, cid_a,
                               sa, afab, sb, dx, bfabs, bct, bcloc, bho,
                               flag_a, vfrac_a, AMREX_D_DECL(apx,apy,apz),
                               AMREX_D_DECL(fcx,fcy,fcz), barea_a, bcent_a, beb);
            });
            BaseFab<GpuArray<Long,eb_stencil>> tmpcolfab
                (bx, 1, (GpuArray<Long,eb_stencil>*)colp);
            amrex::fill(tmpcolfab,
            [=] AMREX_GPU_HOST_DEVICE (GpuArray<Long,eb_stencil>& sten, int i, int j, int k)
            {
                habec_cols_eb(sten, i, j, k, cid_a);
            });
        }
#endif
    }

    Long* last = row_offset.data() + m_nrows_proc;
    ParallelFor(1, [=] AMREX_GPU_DEVICE (Long) noexcept { *last = nentries; });
    Gpu::streamSynchronize();

    if (hasDuplicateColumns(m_geom)) {
        mergeDuplicateColumns(mat, cols, row_offset, m_nrows_proc);
    }

    m_A.define(m_part, mat.data(), cols.data(), nentries, row_offset.data(),
               CsrSorted{false}, CsrValid{false});
}

void
MLAlgMG::Impl::loadRHS (MultiFab const& rhs)
{
    BL_PROFILE("MLAlgMG::loadRHS()");

    if (m_nodal) {
        MLNodeLinOp const* nodelinop = m_nodelinop;
        for (MFIter mfi(m_lid, MFItInfo{}.UseDefaultStream()); mfi.isValid(); ++mfi) {
            if (m_nrows_grid[mfi] == 0) { continue; }
            nodelinop->fillRHS(m_mglev, mfi, m_lid.const_array(mfi),
                               m_b.data() + m_row_begin[mfi], rhs.const_array(mfi));
        }
        Gpu::streamSynchronize();
    } else {
        // Cut-cell rows are weighted by the volume fraction like the matrix
        // rows; zero rhs for overset and covered cells.
        for (MFIter mfi(m_gid, MFItInfo{}.UseDefaultStream()); mfi.isValid(); ++mfi) {
            if (m_nrows_grid[mfi] == 0) { continue; }
            const Box& bx = mfi.validbox();
            Real* bp = m_b.data() + m_row_begin[mfi];
            auto const& rhs_a = rhs.const_array(mfi);
            auto osm = m_overset_mask ? m_overset_mask->const_array(mfi) : Array4<int const>();
#if defined(AMREX_USE_EB) && (AMREX_SPACEDIM > 1)
            bool const cut = m_flags && (*m_flags)[mfi].getType(bx) == FabType::singlevalued;
            auto flag = cut ? m_flags->const_array(mfi) : Array4<EBCellFlag const>();
            auto vfrc = cut ? m_vfrac->const_array(mfi) : Array4<Real const>();
#endif
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                bool norow = (osm && osm(i,j,k) == 0);
                Real w = Real(1.0);
#if defined(AMREX_USE_EB) && (AMREX_SPACEDIM > 1)
                if (flag) {
                    norow = norow || flag(i,j,k).isCovered();
                    w = vfrc(i,j,k);
                }
#endif
                bp[bx.index(IntVect{AMREX_D_DECL(i,j,k)})] = norow ? Real(0.0) : rhs_a(i,j,k) * w;
            });
        }
        Gpu::streamSynchronize();
    }
}

void
MLAlgMG::Impl::getSolution (MultiFab& soln)
{
    BL_PROFILE("MLAlgMG::getSolution()");

    if (m_nodal) {
        soln.setVal(Real(0.0));
        m_tmp.setVal(Real(0.0));
        for (MFIter mfi(m_tmp, MFItInfo{}.UseDefaultStream()); mfi.isValid(); ++mfi) {
            if (m_nrows_grid[mfi] == 0) { continue; }
            const Box& bx = mfi.validbox();
            auto const& xfab = m_tmp.array(mfi);
            auto const& lid = m_lid.const_array(mfi);
            Real const* xp = m_x.data() + m_row_begin[mfi];
            AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
            {
                if (lid(i,j,k) >= 0) {
                    xfab(i,j,k) = xp[lid(i,j,k)];
                }
            });
        }
        // Shared nodes that this rank does not own get the owner's value.
        soln.ParallelAdd(m_tmp, 0, 0, 1, m_geom.periodicity());
    } else {
        soln.setBndry(Real(0.0));
        for (MFIter mfi(m_gid, MFItInfo{}.UseDefaultStream()); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            auto const& s = soln.array(mfi);
            if (m_nrows_grid[mfi] == 0) {
                AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k, { s(i,j,k) = Real(0.0); });
            } else {
                Real const* xp = m_x.data() + m_row_begin[mfi];
                AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
                {
                    s(i,j,k) = xp[bx.index(IntVect{AMREX_D_DECL(i,j,k)})];
                });
            }
        }
        Gpu::streamSynchronize();
    }
}

}
